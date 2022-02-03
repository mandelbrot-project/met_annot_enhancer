# Load packages

import pandas as pd
import numpy as np
import zipfile
import glob
import os
import sys
import time
import shlex
import subprocess
from tqdm import tqdm
from tqdm import tqdm_notebook
from opentree import OT
import json
from pandas import json_normalize
import yaml
import spectral_lib_matcher

# for debug ony shuld be commented later 
from pathlib import Path
p = Path(__file__).parents[2]
print(p)
os.chdir(p)

# Defining display options

pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 100)

# We deactivate the iloc warning see https://stackoverflow.com/a/20627316
pd.options.mode.chained_assignment = None  # default='warn'

# Loading the parameters from yaml file
if not os.path.exists('configs/user_defined/default.yaml'):
    print('No configs/user_defined/default.yaml: copy from configs/default/default.yaml and modifiy according to your needs')
with open (r'configs/user_defined/default.yaml') as file:    
    params_list = yaml.load(file, Loader=yaml.FullLoader)

download_gnps_job = params_list['options'][0]['download_gnps_job']
do_spectral_match = params_list['options'][1]['do_spectral_match']
do_taxo_resolving = params_list['options'][2]['do_taxo_resolving']
keep_lowest_taxon = params_list['options'][3]['keep_lowest_taxon']
output_plots = params_list['options'][4]['output_plots']


gnps_job_id = params_list['paths'][0]['gnps_job_id']
input_folder = params_list['paths'][1]['input_folder']
project_name = params_list['paths'][2]['project_name']
output_folder = params_list['paths'][3]['output_folder']
metadata_path = params_list['paths'][4]['metadata_path']
db_file_path = params_list['paths'][5]['db_file_path']
adducts_pos_path = params_list['paths'][6]['adducts_pos_path']
adducts_neg_path = params_list['paths'][7]['adducts_neg_path']

parent_mz_tol = params_list['spectral_match_params'][0]['parent_mz_tol']
msms_mz_tol = params_list['spectral_match_params'][1]['msms_mz_tol']
min_cos = params_list['spectral_match_params'][2]['min_cos']
min_peaks = params_list['spectral_match_params'][3]['min_peaks']

Top_N_Sample = params_list['repond_params'][0]['Top_N_Sample']
top_to_output= params_list['repond_params'][1]['top_to_output']
ppm_tol = params_list['repond_params'][2]['ppm_tol']
polarity = params_list['repond_params'][3]['polarity']
organism_header = params_list['repond_params'][4]['organism_header']
sampletype_header = params_list['repond_params'][5]['sampletype_header']
use_post_taxo = params_list['repond_params'][6]['use_post_taxo'] # Note we dont use this variable see why !
top_N_chemical_consistency = params_list['repond_params'][7]['top_N_chemical_consistency']
file_extension = params_list['repond_params'][8]['file_extension']
msfile_suffix = params_list['repond_params'][9]['msfile_suffix']
min_score_taxo_ms1 = params_list['repond_params'][10]['min_score_taxo_ms1']
min_score_chemo_ms1 = params_list['repond_params'][11]['min_score_chemo_ms1']

drop_blanks = params_list['plotting_params'][0]['drop_blanks']
multi_plot = params_list['plotting_params'][1]['multi_plot']



path_to_folder = os.path.expanduser(os.path.join(input_folder , gnps_job_id))
path_to_file = os.path.expanduser(os.path.join(input_folder , gnps_job_id + '.zip'))


# Adding expanduser option to expand home path if encoded in the params file
# Fetching the GNPS folder

path_to_gnps_folder = os.path.expanduser(os.path.join(input_folder , gnps_job_id))


if not os.path.exists(path_to_gnps_folder):
    print('No GNPS input folder found : please check the config.yaml and make sure the paths are set correctly')


quantification_table_reformatted_path = os.path.join(path_to_gnps_folder,'quantification_table_reformatted','')


path_to_results_folders =  os.path.expanduser(os.path.join(output_folder, project_name +'/'))
if not os.path.exists(path_to_results_folders):
    os.makedirs(path_to_results_folders)


query_file_path = os.path.join(path_to_gnps_folder,'spectra/specs_ms.mgf')

spectral_match_results_filename = project_name + '_spectral_match_results.tsv'
isdb_results_path = os.path.join(path_to_results_folders, spectral_match_results_filename)

spectral_match_results_repond_filename = project_name + '_spectral_match_results_repond.tsv'
isdb_results_repond_path = os.path.join(path_to_results_folders, spectral_match_results_repond_filename)

spectral_match_results_repond_flat_filename = project_name + '_spectral_match_results_repond_flat.tsv'
isdb_results_repond_flat_path = os.path.join(path_to_results_folders, spectral_match_results_repond_flat_filename)


spectral_match_results_repond_flat_sel_filename = project_name + '_spectral_match_results_repond_flat_selected.tsv'
isdb_results_repond_flat_sel_path = os.path.join(path_to_results_folders, spectral_match_results_repond_flat_sel_filename)





sunburst_chem_filename = project_name + '_chemo_sunburst.html'
sunburst_organisms_filename = project_name + '_organisms_sunburst.html'

treemap_chemo_counted_filename = project_name + '_chemo_treemap_counted.html'
treemap_chemo_intensity_filename = project_name + '_chemo_treemap_intensity.html'

treemap_chemo_multi_counted_filename = project_name + '_chemo_treemap_multi_counted.html'
treemap_chemo_multi_intensity_filename = project_name + '_chemo_treemap_multi_intensity.html'

pivot_table_filename = project_name + '_pivot_table.html'


sunburst_chem_results_path = os.path.join(path_to_results_folders,sunburst_chem_filename)
sunburst_organisms_results_path = os.path.join(path_to_results_folders,sunburst_organisms_filename)

treemap_chemo_counted_results_path = os.path.join(path_to_results_folders,treemap_chemo_counted_filename)
treemap_chemo_intensity_results_path = os.path.join(path_to_results_folders,treemap_chemo_intensity_filename)

treemap_chemo_multi_counted_results_path = os.path.join(path_to_results_folders,treemap_chemo_multi_counted_filename)
treemap_chemo_multi_intensity_results_path = os.path.join(path_to_results_folders,treemap_chemo_multi_intensity_filename)

pivot_table_results_path = os.path.join(path_to_results_folders,pivot_table_filename)



# We load the df flat fromn previously generated annotation results

df4cyto_flat = pd.read_csv(isdb_results_repond_flat_path, sep='\t')


df4cyto_flat_sub = df4cyto_flat.head(150)

# this step is required to avoid the empty leaves error see (https://plotly.com/python/sunburst-charts/#rectangular-data-with-missing-values)
# Strangely enough it was not observed when directly plotting the snburst (maybe because of a different hanflinh of the nan)

df4cyto_flat_sub = df4cyto_flat_sub.replace({np.nan:'None'})


if output_plots == True:

    print('''
    Generating plots... check your web browser !
    ''')

    import plotly.express as px
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go

    fig = px.sunburst(df4cyto_flat_sub, path=['structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class'])
    fig.update_layout(
        #font_family="Courier New",
        title_font_family="Courier New",
        title_font_color="black",
        title_font_size=14,
        legend_title_font_color="black",
        title_text="<b> Overview of the consensus chemical annotions <br> at the NP Classifier pathway, superclass and class level for <br>" + project_name + "</b>",
        title_x=0.5
    )

    fig.update_layout(
        title={
            'text': "<b> Overview of the consensus chemical annotions <br> at the NP Classifier pathway, superclass and class level for <br>" + '<span style="font-size: 20px;">' + project_name + '</span>' + "</b>",
            'y':0.96,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'})

    fig.update_layout(margin=dict(l=50, r=50, t=100, b=50)
    #,paper_bgcolor="Black"
    )

    fig.show()

    fig.write_html(sunburst_chem_results_path,
                full_html=False,
                include_plotlyjs='cdn')

    if keep_lowest_taxon == False :


        fig = px.sunburst(df4cyto_flat, path=['organism_taxonomy_01domain', 'organism_taxonomy_02kingdom', 'organism_taxonomy_03phylum',
                    'organism_taxonomy_04class', 'organism_taxonomy_05order', 'organism_taxonomy_06family', 'organism_taxonomy_07tribe', 'organism_taxonomy_08genus', 'organism_taxonomy_09species', 'organism_taxonomy_10varietas'],
                        )
        fig.update_layout(
            #font_family="Courier New",
            title_font_family="Courier New",
            title_font_color="black",
            title_font_size=14,
            legend_title_font_color="black",
            title_text="<b> Overview of the source organisms of the chemical annotation <br> at the domain, kingdom, phylum, class, order, family, tribe, genus, species and varietas level for <br>" + project_name + "</b>",
            title_x=0.5
        )

        fig.update_layout(
            title={
                'text': "<b> Overview of the source organisms of the chemical annotation <br> at the domain, kingdom, phylum, class, order, family, tribe, genus, species and varietas level for <br>" + '<span style="font-size: 20px;">' + project_name + '</span>' + "</b>",
                'y':0.96,
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top'})

        fig.update_layout(margin=dict(l=50, r=50, t=100, b=50)
        #,paper_bgcolor="Black"
        )

        fig.show()

        fig.write_html(sunburst_organisms_results_path,
                    full_html=False,
                    include_plotlyjs='cdn')


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
'structure_taxonomy_npclassifier_03class_consensus', 'score_input', 'libname',
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


