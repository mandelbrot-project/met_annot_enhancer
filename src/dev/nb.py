# Required libraries
import pandas as pd
import yaml
import json
from pandas import json_normalize
import time


# Loading helpers functions

from helpers import *

# Loading loaders (ahuum) and formatters functions

from loaders import *
from formatters import *

# Loading other functions

import spectral_lib_matcher
from ms1_matcher import ms1_matcher
from taxo_resolver import *
from reponderation_functions import *


# for debug ony should be commented later 
from pathlib import Path
p = Path(__file__).parents[2]
print(p)
os.chdir(p)


# Loading the parameters from yaml file


if not os.path.exists('configs/user_defined/default_nb.yaml'):
    print('No configs/user_defined/default.yaml: copy from configs/default/default.yaml and modify according to your needs')
with open (r'configs/user_defined/default_nb.yaml') as file:    
    params_list = yaml.load(file, Loader=yaml.FullLoader)

# Parameters can now be accessed using params_list['level1']['level2'] e. g. arams_list['options']['download_gnps_job']

# Downloading GNPS files
if params_list['options']['download_gnps_job'] == True:

    gnps_job_fetcher(gnps_job_id = params_list['paths']['gnps_job_id'], input_folder = params_list['paths']['input_folder'])

# Generating pathes
# The pathes are stored in a dictionary and can then be accesed by paths_dic['value']

paths_dic = paths_generator(params_list = params_list)



# Writing used parameters 
params_suffix = '.yaml'

with open(os.path.join(paths_dic['path_to_results_folders'], params_list['paths']['gnps_job_id'] + params_suffix), 'w') as file:  
    documents = yaml.dump(params_list, file)

print('''
Parameters used are stored in '''
+ str(os.path.join(paths_dic['path_to_results_folders'], params_list['paths']['gnps_job_id'] + params_suffix))
)


# timer is started
start_time = time.time()


######################################################################################################
######################################################################################################
# Spectral matching stage
# If required the msms spectral matching is done 
######################################################################################################


if params_list['options']['do_spectral_match'] == True:

    print('''
    Proceeding to spectral matching ...
    ''')

    spectral_lib_matcher.main(query_file_path=paths_dic['query_file_path'],
                              db_file_path=params_list['paths']['db_file_path'],
                              parent_mz_tol=params_list['spectral_match_params']['parent_mz_tol'],
                              msms_mz_tol=params_list['spectral_match_params']['msms_mz_tol'],
                              min_cos=params_list['spectral_match_params']['min_cos'],
                              min_peaks=params_list['spectral_match_params']['min_peaks'],
                              output_file_path=paths_dic['isdb_results_path']
                              )

    print('''
    Spectral matching done !
    ''')

######################################################################################################


######################################################################################################
######################################################################################################
# Loading of external tables
######################################################################################################
# here we will load all required tables

isdb_results = isdb_results_loader(
    isdb_results_path=paths_dic['isdb_results_path'])

clusterinfo_summary = clusterinfo_summary_loader(
    clusterinfo_summary_path=paths_dic['clusterinfo_summary_path'])

dt_isdb_metadata = isdb_metadata_loader(
    isdb_metadata_path=params_list['paths']['metadata_path'])

dt_samples_metadata = samples_metadata_loader(samples_metadata_table_path=paths_dic['samples_metadata_table_path'],
                                              organism_header=params_list['repond_params']['organism_header'])

feature_intensity_table = feature_intensity_table_loader(
    feature_intensity_table_path=paths_dic['quantification_table_reformatted_path'])

######################################################################################################
######################################################################################################
# Merging
######################################################################################################
# We want to complement the ISDB results file with the component index and 
# the parent mass of the ion (for this last one this could be done earlier)


dt_isdb_results = pd.merge(isdb_results, clusterinfo_summary, on='feature_id')


print('Number of features: ' + str(len(clusterinfo_summary)))
print('Number of MS2 annotation: ' + str(len(dt_isdb_results)))


######################################################################################################
######################################################################################################
# MS1 matching stage
######################################################################################################
# Now we directly do the MS1 matching stage on the cluster_summary. 
# No need to have MS2 annotations


df_MS1_matched = ms1_matcher(input_df=clusterinfo_summary,
                             adducts_file_path=params_list['paths']['adducts_pos_path'],
                             ppm_tol=params_list['repond_params']['ppm_tol'],
                             df_metadata=dt_isdb_metadata)



# ######################################################################################################
# ######################################################################################################
# # Merging
# ######################################################################################################
# # We can now merge MS1 and MS2 annotations


dt_isdb_results = pd.concat([dt_isdb_results, df_MS1_matched])


print('Number of annotated features at the MS1 level : ' +
      str(len(df_MS1_matched['feature_id'].unique())))

print('Total number of unique MS1 and MSMS annotations: ' + str(len(dt_isdb_results)))




######################################################################################################
######################################################################################################
# Merging
######################################################################################################
# We now complement the previous list of annotation with selected fields of the metadata table


# # Rank annotations based on the spectral score

dt_isdb_results["msms_score"] = pd.to_numeric(
    dt_isdb_results["msms_score"], downcast="float")
dt_isdb_results['rank_spec'] = dt_isdb_results[['feature_id', 'msms_score']].groupby(
    'feature_id')['msms_score'].rank(method='dense', ascending=False)

dt_isdb_results.reset_index(inplace=True, drop=True)

# now we merge with the Occurences DB metadata after selection of our columns of interest

cols_to_use = ['structure_inchikey', 'structure_inchi',
            'structure_smiles', 'structure_molecular_formula',
            'structure_exact_mass', 'short_inchikey', 'structure_taxonomy_npclassifier_01pathway', 
            'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class',
            'organism_name', 'organism_taxonomy_ottid',
            'organism_taxonomy_01domain', 'organism_taxonomy_02kingdom', 'organism_taxonomy_03phylum',
            'organism_taxonomy_04class', 'organism_taxonomy_05order', 'organism_taxonomy_06family', 'organism_taxonomy_07tribe', 'organism_taxonomy_08genus', 'organism_taxonomy_09species', 'organism_taxonomy_10varietas' ]

dt_isdb_results = pd.merge(
    left=dt_isdb_results, right=dt_isdb_metadata[cols_to_use], left_on='short_inchikey', right_on='short_inchikey', how='outer')
dt_isdb_results.dropna(subset=['feature_id'], inplace=True)

print('Total number of annotations with unique biosource per line: ' +
      str(len(dt_isdb_results)))



######################################################################################################
######################################################################################################
# Taxonomical resolving
######################################################################################################
# Resolving the taxon information from the samples metadata file


dt_samples_metadata = taxa_lineage_appender(samples_metadata=dt_samples_metadata,
                                         organism_header=params_list['repond_params']['organism_header'],
                                         do_taxo_resolving=params_list['options']['do_taxo_resolving'],
                                         path_to_results_folders=paths_dic['path_to_results_folders'],
                                         project_name=params_list['paths']['project_name'])




######################################################################################################
######################################################################################################
# Establishing the contribution of individual biosources for each feature
######################################################################################################

# First the feature table intensity file is loaded and appropriately formatted

feature_intensity_table_formatted = feature_intensity_table_formatter(feature_intensity_table_path=paths_dic['quantification_table_reformatted_path'],
                                                                      file_extension=params_list[
                                                                          'repond_params']['file_extension'],
                                                                      msfile_suffix=params_list['repond_params']['msfile_suffix'])


# The function below will fetch the topN most contributing biosources for each feature and outputs them as a dataframe of lists 


topN_contributors = biosource_contribution_fetcher(feature_intensity_table= feature_intensity_table_formatted,
                               samples_metadata= dt_samples_metadata ,
                               top_n=params_list['repond_params']['Top_N_Sample'])


# The top N contributors list is now appended to the annotation table 

dt_isdb_results_topN = pd.merge(
    dt_isdb_results, topN_contributors, left_on='feature_id', right_on='row_ID', how='left')


######################################################################################################
######################################################################################################
# Taxonomical Reweighting
######################################################################################################

dt_taxo_reweighed = taxonomical_reponderator(dt_isdb_results=dt_isdb_results_topN,
                                             min_score_taxo_ms1=params_list['repond_params']['min_score_taxo_ms1'])


######################################################################################################
######################################################################################################
# Structural Consistency Reweighting
######################################################################################################


dt_taxo_chemo_reweighed = chemical_reponderator(clusterinfo_summary_file=clusterinfo_summary,
                                                dt_isdb_results=dt_taxo_reweighed,
                                                top_N_chemical_consistency=params_list['repond_params']['top_N_chemical_consistency'])


######################################################################################################
######################################################################################################
# Filtering
######################################################################################################


dt_taxo_chemo_reweighed_topN = top_N_slicer(dt_isdb_results=dt_taxo_chemo_reweighed,
                                            top_to_output=params_list['repond_params']['top_to_output'])


dt_taxo_chemo_reweighed_topN_filt = annotation_table_formatter(dt_input=dt_taxo_chemo_reweighed_topN,
                                                               keep_lowest_taxon=params_list['options']['keep_lowest_taxon'],
                                                               min_score_taxo_ms1=params_list[
                                                                   'repond_params']['min_score_taxo_ms1'],
                                                               min_score_chemo_ms1=params_list['repond_params']['min_score_chemo_ms1'])


######################################################################################################
######################################################################################################
# Fetching CHEMBL Ids
######################################################################################################

dt_taxo_chemo_reweighed_topN_filt = chembl_id_fetcher(df_input = dt_taxo_chemo_reweighed_topN_filt )


all_columns = list(df4cyto_flat) # Creates list of all column headers
df4cyto_flat[all_columns] = df4cyto_flat[all_columns].astype(str)

gb_spec = {c: '|'.join for c in annot_attr}
for c in comp_attr:
    gb_spec[c] = 'first'

df4cyto = df4cyto_flat.groupby('feature_id').agg(gb_spec)








df4cyto_flat.to_csv(isdb_results_repond_flat_path, sep='\t')

df4cyto.to_csv(isdb_results_repond_path, sep='\t')





print('Finished in %s seconds.' % (time.time() - start_time))
print('You can check your results here %s' % isdb_results_repond_path)





if output_plots == True:

    print('''
    Generating plots... check your web browser !
    ''')

    import plotly.express as px
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go


    # if keep_lowest_taxon == False :


    #     # we have a problem because the organism_taxonomy_ are lists and not strings.
    #     # We subset specifically these columns

    #     dt_isdb_results_tax = dt_isdb_results.loc[:, dt_isdb_results.columns.str.startswith('organism_taxonomy_')]

    #     # and then use the explode function to yield the datframe with the values extractes from the lists

    #     dt_isdb_results_tax = dt_isdb_results_tax.set_index(['organism_taxonomy_ottid']).apply(pd.Series.explode).reset_index()

    #     # we now drop the previous columns with the list format 

    #     colsToDrop = [ 'organism_taxonomy_01domain', 'organism_taxonomy_02kingdom',
    #         'organism_taxonomy_03phylum', 'organism_taxonomy_04class', 'organism_taxonomy_05order',
    #         'organism_taxonomy_06family', 'organism_taxonomy_07tribe', 'organism_taxonomy_08genus',
    #         'organism_taxonomy_09species', 'organism_taxonomy_10varietas']

    #     dt_isdb_results = dt_isdb_results.drop(colsToDrop, axis=1)

    #     # and we merge back using the organism_taxonomy_ottid. Here we use concat and merge on indexes

    #     dt_isdb_results = pd.concat([dt_isdb_results, dt_isdb_results_tax], axis=1)



    #     dt_isdb_results['counter'] = 1

    #     dt_isdb_results = dt_isdb_results.replace({np.nan:'None'})

    #     fig = px.treemap(dt_isdb_results, path=[px.Constant("all"), 'organism_taxonomy_01domain', 'organism_taxonomy_02kingdom', 'organism_taxonomy_03phylum',
    #                 'organism_taxonomy_04class', 'organism_taxonomy_05order', 'organism_taxonomy_06family', 'organism_taxonomy_07tribe', 'organism_taxonomy_08genus', 'organism_taxonomy_09species'],  values='counter')


    #     fig.show()

    #     fig.update_layout(
    #         title_font_family="Courier New",
    #         title_font_color="black",
    #         title_font_size=14,
    #         legend_title_font_color="black",
    #         title_text="<b> Overview of the source organisms taxonomical repartition of the chemical annotations <br> before taxonomical reponderation <br>" + project_name + "</b>",
    #         title_x=0.5
    #     )

    #     fig.update_layout(
    #         title={
    #             'text': "<b> Overview of the source organisms taxonomical repartition of the chemical annotations <br> before taxonomical reponderation for <br>" + '<span style="font-size: 20px;">' + project_name + '</span>' + "</b>",
    #             'y':0.96,
    #             'x':0.5,
    #             'xanchor': 'center',
    #             'yanchor': 'top'})

    #     fig.update_layout(margin=dict(l=50, r=50, t=100, b=50)
    #     #,paper_bgcolor="Black"
    #     )

    #     fig.show()

        # fig.write_html(sunburst_organisms_results_path,
        #             full_html=False,
        #             include_plotlyjs='cdn')


    # here we want to have an puput per sample so we merge back the annotation frame with the feature table 
    #actually we might want to have the metadat joined to the feature table since the beginning

    feature_intensity_table_t

    samples_metadata_full = pd.read_csv(metadata_table_path + str(os.listdir(metadata_table_path)[0]), sep='\t')
    feature_intensity_meta = pd.merge(left=samples_metadata_full, right=feature_intensity_table_t, left_on='filename', right_on='MS_filename',how='inner')


    feature_intensity_meta_gp_species = feature_intensity_meta.groupby(organism_header).mean()
    feature_intensity_meta_gp_species = feature_intensity_meta_gp_species.transpose()
    feature_intensity_meta_gp_species.index.name = 'row_ID'



    feature_intensity_table = feature_intensity_table_t.transpose()

    feature_intensity_table.reset_index(inplace=True)
    feature_intensity_meta_gp_species.reset_index(inplace=True)

    

    ft_merged = pd.merge(feature_intensity_table, feature_intensity_meta_gp_species, on='row_ID', how='left')

    if multi_plot == True:
        feature_intensity_meta_gp_multi = feature_intensity_meta.groupby([organism_header,sampletype_header]).mean()
        feature_intensity_meta_gp_multi = feature_intensity_meta_gp_multi.transpose()
        feature_intensity_meta_gp_multi.columns = feature_intensity_meta_gp_multi.columns.map('_'.join)
        feature_intensity_meta_gp_multi.index.name = 'row_ID'
        feature_intensity_meta_gp_multi.reset_index(inplace=True)
        ft_merged = pd.merge(ft_merged, feature_intensity_meta_gp_multi, on='row_ID', how='left')

    df4cyto_flat['feature_id'] = df4cyto_flat['feature_id'].astype('int')

    dt_isdb_results_int = pd.merge(
        df4cyto_flat, ft_merged, left_on='feature_id', right_on='row_ID', how='left')

    dt_isdb_results_int['counter'] = 1




    # for n in samples_metadata_full['species_name'].unique():
    #     print(n)

    #     dt = dt_isdb_results_int[dt_isdb_results_int[n] > 0]

    #     fig = px.treemap(dt, path=[px.Constant("all"), 'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class'], values=n)
    #     fig.update_traces(root_color="lightgrey")
    #     fig.update_layout(margin = dict(t=50, l=25, r=25, b=25))
    #     fig.show()


    # for n in samples_metadata_full['species_name'].unique():
    #     print(n)

    #     dt = dt_isdb_results_int[dt_isdb_results_int[n] > 0]

    #     fig = px.treemap(dt, path=[px.Constant("all"), 'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class'], values='counter')
    #     fig.update_traces(root_color="lightgrey")
    #     fig.update_layout(margin = dict(t=50, l=25, r=25, b=25))
    #     fig.show()



    # here we also need to incresae the spec list of list of dic according to the lenght of unique_group_labels this can be done following https://stackoverflow.com/a/3459131 
    # in fact we do this with itertools
    if drop_blanks == True:
        samples_metadata_full = samples_metadata_full[~samples_metadata_full[sampletype_header].str.contains("none|BK|blanck|bk|mock")]

    import itertools
    unique_group_labels = samples_metadata_full[organism_header].unique()

    pattern=[{"type": "domain"}]


    rep_pattern = list(itertools.chain.from_iterable(itertools.repeat(x, len(unique_group_labels)) for x in pattern))


    fig = make_subplots(1, len(unique_group_labels),
    subplot_titles = (unique_group_labels),
    specs=[rep_pattern])

    i=1
    for n in unique_group_labels:
        print(n)

        dt = dt_isdb_results_int[dt_isdb_results_int[n] > 0]
        fig.add_trace(px.treemap(dt, path=[px.Constant("all"), 'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class'], 
        values='counter').data[0], 
        row=1,col=i)
        i+=1

    fig.update_traces(root_color="lightgrey")
    fig.update_layout(margin = dict(t=50, l=25, r=25, b=25),
    title_text="Metabolite annotation overview (size proportional to individual count)")
    fig.update_annotations(font_size=12)
    fig.show()
    fig.write_html(treemap_chemo_counted_results_path,
                full_html=False,
                include_plotlyjs='cdn')

    fig = make_subplots(1, len(unique_group_labels),
    subplot_titles = (unique_group_labels),
    specs=[rep_pattern])

    i=1
    for n in unique_group_labels:
        print(n)

        dt = dt_isdb_results_int[dt_isdb_results_int[n] > 0]
        fig.add_trace(px.treemap(dt, path=[px.Constant("all"), 'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class'], 
        values=n).data[0], 
        row=1,col=i)
        i+=1

    fig.update_traces(root_color="lightgrey")
    fig.update_layout(margin = dict(t=50, l=25, r=25, b=25),
    title_text="Metabolite annotation overview (size proportional to mean intensity)")
    fig.update_annotations(font_size=12)
    fig.show()
    fig.write_html(treemap_chemo_intensity_results_path,
                full_html=False,
                include_plotlyjs='cdn')

#### working on the multilabelled subset

    if multi_plot == True:
        
        samples_metadata_full['combined'] = samples_metadata_full[organism_header] + '_' + samples_metadata_full[sampletype_header]
        unique_group_labels = samples_metadata_full['combined'].unique()
        type(unique_group_labels)

        pattern=[{"type": "domain"}]


        rep_pattern = list(itertools.chain.from_iterable(itertools.repeat(x, len(unique_group_labels)) for x in pattern))

        fig = make_subplots(1, len(unique_group_labels),
        subplot_titles = (unique_group_labels),
        specs=[rep_pattern])

        i=1
        for n in unique_group_labels:
            print(n)

            dt = dt_isdb_results_int[dt_isdb_results_int[n] > 0]
            fig.add_trace(px.treemap(dt, path=[px.Constant("all"), 'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class'], 
            values='counter').data[0], 
            row=1,col=i)
            i+=1

        fig.update_traces(root_color="lightgrey")
        fig.update_layout(margin = dict(t=50, l=25, r=25, b=25),
        title_text="Metabolite annotation overview (size proportional to individual count)")
        fig.update_annotations(font_size=12)
        fig.show()
        fig.write_html(treemap_chemo_multi_counted_results_path,
                    full_html=False,
                    include_plotlyjs='cdn')



        fig = make_subplots(1, len(unique_group_labels),
        subplot_titles = (unique_group_labels),
        specs=[rep_pattern])

        i=1
        for n in unique_group_labels:
            print(n)

            dt = dt_isdb_results_int[dt_isdb_results_int[n] > 0]
            fig.add_trace(px.treemap(dt, path=[px.Constant("all"), 'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class'], 
            values=n).data[0], 
            row=1,col=i)
            i+=1

        fig.update_traces(root_color="lightgrey")
        fig.update_layout(margin = dict(t=50, l=25, r=25, b=25),
        title_text="Metabolite annotation overview (size proportional to mean intensity)")
        fig.update_annotations(font_size=12)
        fig.show()
        fig.write_html(treemap_chemo_multi_intensity_results_path,
                    full_html=False,
                    include_plotlyjs='cdn')


df4cyto_flat['final_score'] = df4cyto_flat['final_score'].astype('float')
df4cyto_flat[df4cyto_flat['final_score'] >= 8]

df4cyto_flat.columns


df4cyto_flat_sel = df4cyto_flat[['feature_id', 'component_id', 'structure_taxonomy_npclassifier_01pathway_consensus','structure_taxonomy_npclassifier_02superclass_consensus',
'structure_taxonomy_npclassifier_03class_consensus', 'msms_score', 'libname',
'structure_inchikey', 'structure_inchi', 'structure_smiles', 'structure_molecular_formula',
'adduct', 'structure_exact_mass', 'short_inchikey',
'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass',
'structure_taxonomy_npclassifier_03class', 'organism_name', 'organism_taxonomy_ottid',
'organism_taxonomy_01domain', 'organism_taxonomy_02kingdom', 'organism_taxonomy_03phylum',
'organism_taxonomy_04class', 'organism_taxonomy_05order', 'organism_taxonomy_06family',
'organism_taxonomy_07tribe', 'organism_taxonomy_08genus', 'organism_taxonomy_09species',
'score_taxo', 'score_max_consistency', 'final_score']]



df4cyto_flat_sel.to_csv(isdb_results_repond_flat_sel_path, sep='\t', index=None)


from pivottablejs import pivot_ui

pivot_ui(dt_isdb_results_int, outfile_path=pivot_table_results_path)


