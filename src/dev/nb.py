# %%
# Required libraries
import pandas as pd
import yaml
import json
from pandas import json_normalize
import time
# %%
# Loading helpers functions

from helpers import *
# %%
# Loading loaders (ahuum) and formatters functions

from loaders import *
from formatters import *
# %%
# Loading other functions

import spectral_lib_matcher
from ms1_matcher import ms1_matcher
from taxo_resolver import *
from reponderation_functions import *
from plotter import *

# %%
# for debug ony should be commented later
from pathlib import Path
p = Path(__file__).parents[2]
print(p)
os.chdir(p)

# %%
# Loading the parameters from yaml file


if not os.path.exists('configs/user_defined/default.yaml'):
    print('No configs/user_defined/default.yaml: copy from configs/default/default.yaml and modify according to your needs')
with open(r'configs/user_defined/default.yaml') as file:
    params_list = yaml.load(file, Loader=yaml.FullLoader)

# Parameters can now be accessed using params_list['level1']['level2'] e. g. arams_list['options']['download_gnps_job']

# %%
# Downloading GNPS files
if params_list['options']['download_gnps_job'] == True:

    gnps_job_fetcher(gnps_job_id=params_list['paths']['gnps_job_id'],
                     input_folder=params_list['paths']['input_folder'])
# %%
# Generating pathes
# The pathes are stored in a dictionary and can then be accesed by paths_dic['value']

paths_dic = paths_generator(params_list=params_list)


# Writing used parameters
params_suffix = '.yaml'


with open(os.path.join(paths_dic['path_to_results_folders'], params_list['paths']['gnps_job_id'] + params_suffix), 'w') as file:
    documents = yaml.dump(params_list, file)

print('''
Parameters used are stored in '''
      + str(os.path.join(paths_dic['path_to_results_folders'],
            params_list['paths']['gnps_job_id'] + params_suffix))
      )

# %%
# timer is started
start_time = time.time()

# %%
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

# %%
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
    isdb_metadata_path=params_list['paths']['metadata_path'],
    organism_header=params_list['metadata_params']['organism_header'])

dt_samples_metadata = samples_metadata_loader(samples_metadata_table_path=paths_dic['samples_metadata_table_path'],
                                              organism_header=params_list['repond_params']['organism_header'])

feature_intensity_table = feature_intensity_table_loader(
    feature_intensity_table_path=paths_dic['quantification_table_reformatted_path'])

# %%
######################################################################################################
######################################################################################################
# Merging
######################################################################################################
# We want to complement the ISDB results file with the component index and
# the parent mass of the ion (for this last one this could be done earlier)


dt_isdb_results = pd.merge(isdb_results, clusterinfo_summary, on='feature_id')


print('Number of features: ' + str(len(clusterinfo_summary)))
print('Number of MS2 annotation: ' + str(len(dt_isdb_results)))


# %%
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

# %%
# ######################################################################################################
# ######################################################################################################
# # Merging
# ######################################################################################################
# # We can now merge MS1 and MS2 annotations


dt_isdb_results = pd.concat([dt_isdb_results, df_MS1_matched])


print('Number of annotated features at the MS1 level : ' +
      str(len(df_MS1_matched['feature_id'].unique())))

print('Total number of unique MS1 and MSMS annotations: ' +
      str(len(dt_isdb_results)))

# %%
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
               'organism_taxonomy_04class', 'organism_taxonomy_05order', 'organism_taxonomy_06family', 'organism_taxonomy_07tribe', 'organism_taxonomy_08genus', 'organism_taxonomy_09species', 'organism_taxonomy_10varietas']

dt_isdb_results = pd.merge(
    left=dt_isdb_results, right=dt_isdb_metadata[cols_to_use], left_on='short_inchikey', right_on='short_inchikey', how='outer')
dt_isdb_results.dropna(subset=['feature_id'], inplace=True)

print('Total number of annotations with unique biosource per line: ' +
      str(len(dt_isdb_results)))

# %%
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

# %%
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


topN_contributors = biosource_contribution_fetcher(feature_intensity_table=feature_intensity_table_formatted,
                                                   samples_metadata=dt_samples_metadata,
                                                   top_n=params_list['repond_params']['Top_N_Sample'])


# The top N contributors list is now appended to the annotation table

dt_isdb_results_topN = pd.merge(
    dt_isdb_results, topN_contributors, left_on='feature_id', right_on='row_ID', how='left')

# %%
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
# Filtering top N hits
######################################################################################################


dt_taxo_chemo_reweighed_topN = top_N_slicer(dt_isdb_results=dt_taxo_chemo_reweighed,
                                            top_to_output=params_list['repond_params']['top_to_output'])

# %%
######################################################################################################
######################################################################################################
# Fetching CHEMBL Ids
######################################################################################################

dt_taxo_chemo_reweighed_chembl = chembl_id_fetcher(
    df_input=dt_taxo_chemo_reweighed_topN)


######################################################################################################
######################################################################################################
# Final formattings
######################################################################################################


df_flat, df_for_cyto = annotation_table_formatter(dt_input=dt_taxo_chemo_reweighed_chembl,
                                                  keep_lowest_taxon=params_list['options']['keep_lowest_taxon'],
                                                  min_score_taxo_ms1=params_list[
                                                      'repond_params']['min_score_taxo_ms1'],
                                                  min_score_chemo_ms1=params_list['repond_params']['min_score_chemo_ms1'])

# %%
######################################################################################################
######################################################################################################
# Exporting annotation results
######################################################################################################


df_flat.to_csv(paths_dic['isdb_results_repond_flat_path'], sep='\t')

df_for_cyto.to_csv(paths_dic['isdb_results_repond_path'], sep='\t')

# %%
######################################################################################################
######################################################################################################
# Preparing tables for plots
######################################################################################################
# Loading clean tables


feature_intensity_table_formatted = feature_intensity_table_formatter(feature_intensity_table_path=paths_dic['quantification_table_reformatted_path'],
                                                                      file_extension=params_list[
                                                                          'repond_params']['file_extension'],
                                                                      msfile_suffix=params_list['repond_params']['msfile_suffix'])

dt_samples_metadata = samples_metadata_full_loader(
    samples_metadata_table_path=paths_dic['samples_metadata_table_path'])

# Formatting the tables

table_for_plots_formatted = table_for_plots_formatter(df_flat=df_flat,
                                                      feature_intensity_table_formatted=feature_intensity_table_formatted,
                                                      dt_samples_metadata=dt_samples_metadata,
                                                      organism_header=params_list['repond_params']['organism_header'],
                                                      sampletype_header=params_list['repond_params']['sampletype_header'],
                                                      multi_plot=params_list['plotting_params']['multi_plot'])

# Some optional filtering can be done

samples_metadata_filtered = samples_metadata_filterer(dt_samples_metadata=dt_samples_metadata,
                                                      organism_header=params_list['repond_params']['organism_header'],
                                                      sampletype_header=params_list['repond_params']['sampletype_header'],
                                                      drop_pattern=params_list['plotting_params']['drop_pattern'])

# %%
######################################################################################################
######################################################################################################
# Plotting figures
######################################################################################################
# Single parameters


plotter_single(dt_isdb_results_int=table_for_plots_formatted,
               dt_samples_metadata=samples_metadata_filtered,
               organism_header=params_list['repond_params']['organism_header'],
               treemap_chemo_counted_results_path=paths_dic['treemap_chemo_counted_results_path'],
               treemap_chemo_intensity_results_path=paths_dic['treemap_chemo_intensity_results_path'])

if params_list['plotting_params']['multi_plot'] == True:

    plotter_multi(dt_isdb_results_int=table_for_plots_formatted,
                dt_samples_metadata=samples_metadata_filtered,
                organism_header=params_list['repond_params']['organism_header'],
                sampletype_header=params_list['repond_params']['sampletype_header'],
                treemap_chemo_multi_counted_results_path=paths_dic['treemap_chemo_multi_counted_results_path'],
                treemap_chemo_multi_intensity_results_path=paths_dic['treemap_chemo_multi_intensity_results_path'])


pivot_tabler(df_input= dt_taxo_chemo_reweighed_chembl,
             lib_to_keep=params_list['filtering_params']['lib_to_keep'],
             minimal_taxo_score=params_list['filtering_params']['minimal_taxo_score'],
             minimal_chemo_score=params_list['filtering_params']['minimal_chemo_score'],
             minimal_total_score=params_list['filtering_params']['minimal_total_score'],
             isdb_results_repond_flat_sel_path=paths_dic['isdb_results_repond_flat_sel_path'],
             pivot_table_results_path=paths_dic['pivot_table_results_path'])


