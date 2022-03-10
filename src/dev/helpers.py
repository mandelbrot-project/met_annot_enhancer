# required libraries

import pandas as pd
import numpy as np
import zipfile
import glob
import os
import sys
import shlex
import subprocess
import json
from pandas import json_normalize
from chembl_webresource_client.new_client import new_client
from tqdm import tqdm
from pivottablejs import pivot_ui



def gnps_job_fetcher(gnps_job_id, input_folder):
    """Fetch a GNPS job and saves its to a given folder
    Args:
            gnps_job_id (str): a GNPS Job id
            input_folder (str): the folder where the job should be kept
    Returns:
        nothing
    """
    path_to_folder = os.path.expanduser(
        os.path.join(input_folder, gnps_job_id))
    path_to_file = os.path.expanduser(
        os.path.join(input_folder, gnps_job_id + '.zip'))

    print('''
    Fetching the GNPS job: '''
    + gnps_job_id
    )
    # or &view=download_clustered_spectra or download_cytoscape_data (check when to use which and optionalize)
    job_url_zip = "https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task=" + \
        gnps_job_id+"&view=download_cytoscape_data"

    cmd = 'curl -d "" '+job_url_zip+' -o '+path_to_file + ' --create-dirs'
    subprocess.call(shlex.split(cmd))

    with zipfile.ZipFile(path_to_file, 'r') as zip_ref:
        zip_ref.extractall(path_to_folder)

    # We finally remove the zip file
    os.remove(path_to_file)

    print('''
    Job successfully downloaded: results are in: '''
    + path_to_folder
    )


def paths_generator(params_list):


    """Generates pathes used by the script according to parameters of the yaml file
    Args:
        params_list (list) : the yaml parameters
    Returns:
        pathes (str)
    """
    # We set an empty dic which will keep the pathes and their key

    paths_dic = {}
    
    # Adding expanduser option to expand home path if encoded in the params file
    # Fetching the GNPS folder
    input_folder = params_list['paths']['input_folder']
    gnps_job_id = params_list['paths']['gnps_job_id']

    path_to_gnps_folder = os.path.expanduser(os.path.join(input_folder , gnps_job_id))


    if not os.path.exists(path_to_gnps_folder):

        print('No GNPS input folder found : please check the config.yaml and make sure the paths are set correctly')



    output_folder = params_list['paths']['output_folder']
    project_name = params_list['paths']['project_name']

    paths_dic['path_to_results_folders'] =  os.path.expanduser(os.path.join(output_folder, project_name +'/'))
    if not os.path.exists(paths_dic['path_to_results_folders']):
        os.makedirs(paths_dic['path_to_results_folders'])
    
    # Generating GNPS pathes
    paths_dic['clusterinfo_summary_path'] = os.path.join(path_to_gnps_folder,'clusterinfo_summary','')
    paths_dic['samples_metadata_table_path'] = os.path.join(path_to_gnps_folder,'metadata_table','')
    paths_dic['quantification_table_reformatted_path'] = os.path.join(path_to_gnps_folder,'quantification_table_reformatted','')


    paths_dic['query_file_path'] = os.path.join(path_to_gnps_folder,'spectra/specs_ms.mgf')

    spectral_match_results_filename = project_name + '_spectral_match_results.tsv'
    paths_dic['isdb_results_path'] = os.path.join(paths_dic['path_to_results_folders'], spectral_match_results_filename)

    spectral_match_results_repond_filename = project_name + '_spectral_match_results_repond.tsv'
    paths_dic['isdb_results_repond_path'] = os.path.join(paths_dic['path_to_results_folders'], spectral_match_results_repond_filename)

    spectral_match_results_repond_flat_filename = project_name + '_spectral_match_results_repond_flat.tsv'
    paths_dic['isdb_results_repond_flat_path'] = os.path.join(paths_dic['path_to_results_folders'], spectral_match_results_repond_flat_filename)


    spectral_match_results_repond_flat_sel_filename = project_name + '_spectral_match_results_repond_flat_selected.tsv'
    paths_dic['isdb_results_repond_flat_sel_path'] = os.path.join(paths_dic['path_to_results_folders'], spectral_match_results_repond_flat_sel_filename)


    sunburst_chem_filename = project_name + '_chemo_sunburst.html'
    sunburst_organisms_filename = project_name + '_organisms_sunburst.html'

    treemap_chemo_counted_filename = project_name + '_chemo_treemap_counted.html'
    treemap_chemo_intensity_filename = project_name + '_chemo_treemap_intensity.html'

    treemap_chemo_multi_counted_filename = project_name + '_chemo_treemap_multi_counted.html'
    treemap_chemo_multi_intensity_filename = project_name + '_chemo_treemap_multi_intensity.html'

    pivot_table_filename = project_name + '_pivot_table.html'

    paths_dic['sunburst_chem_results_path'] = os.path.join(paths_dic['path_to_results_folders'],sunburst_chem_filename)
    paths_dic['sunburst_organisms_results_path'] = os.path.join(paths_dic['path_to_results_folders'],sunburst_organisms_filename)

    paths_dic['treemap_chemo_counted_results_path'] = os.path.join(paths_dic['path_to_results_folders'],treemap_chemo_counted_filename)
    paths_dic['treemap_chemo_intensity_results_path'] = os.path.join(paths_dic['path_to_results_folders'],treemap_chemo_intensity_filename)

    paths_dic['treemap_chemo_multi_counted_results_path'] = os.path.join(paths_dic['path_to_results_folders'],treemap_chemo_multi_counted_filename)
    paths_dic['treemap_chemo_multi_intensity_results_path'] = os.path.join(paths_dic['path_to_results_folders'],treemap_chemo_multi_intensity_filename)

    paths_dic['pivot_table_results_path'] = os.path.join(paths_dic['path_to_results_folders'],pivot_table_filename)

    return paths_dic

def cluster_counter(clusterinfo_summary_file):


    """ Count the numbers of nodes per component index in a molecular network

    Args:
        clusterinfo_summary_file (dataframe) : a molecular network clusterinfo_summary_file
    Returns:
        cluster_count (dataframe): a dataframe with the number of nodes per component index
    """


    cluster_count = clusterinfo_summary_file.drop_duplicates(
        subset=['feature_id', 'component_id']).groupby("component_id").count()
    cluster_count = cluster_count[['feature_id']].rename(
        columns={'feature_id': 'ci_count'}).reset_index()
    return cluster_count


def top_N_slicer(dt_isdb_results, top_to_output):

    """ Keeps only the top N candidates out of an annotation table and sorts them by rank

    Args:
        dt_isdb_results (dataframe) : a annotation table
        top_to_output (integer): Top N of candidate to keep
    Returns:
        dt_isdb_results_chem_rew (dataframe): a dataframe with the top N annotation ordered by final rank
    """

    dt_isdb_results_chem_rew = dt_isdb_results.loc[(
        dt_isdb_results.rank_final <= int(top_to_output))]
    dt_isdb_results_chem_rew[["feature_id", "rank_final", "component_id"]] = dt_isdb_results_chem_rew[[
        "feature_id", "rank_final", "component_id"]].apply(pd.to_numeric, downcast='signed', axis=1)
    dt_isdb_results_chem_rew = dt_isdb_results_chem_rew.sort_values(
        ["feature_id", "rank_final"], ascending=(False, True))

    return dt_isdb_results_chem_rew


def annotation_table_formatter(dt_input, keep_lowest_taxon, min_score_taxo_ms1, min_score_chemo_ms1):
    """ A bunch of filtering functions

    Args:
        dt_isdb_results (dataframe) : a annotation table
        top_to_output (integer): Top N of candidate to keep
    Returns:
        dt_isdb_results_chem_rew (dataframe): a dataframe with the top N annotation ordered by final rank
    """

    # Here we would like to filter results when short IK are repeated for the same feature_id at the same final rank

    dt_input = dt_input.drop_duplicates(
        subset=['feature_id', 'short_inchikey'], keep='first')

    dt_input = dt_input.astype(
        {'feature_id': 'int64'})

    if keep_lowest_taxon == True:

        dt_input['lowest_matched_taxon'] = dt_input['matched_species']
        dt_input['lowest_matched_taxon'] = dt_input['lowest_matched_taxon'].replace(
            'nan', np.NaN)
        col_matched = ['matched_genus', 'matched_tribe', 'matched_family', 'matched_order',
                       'matched_order', 'matched_phylum', 'matched_kingdom', 'matched_domain']
        for col in col_matched:
            dt_input[col] = dt_input[col].replace(
                'nan', np.NaN)
            dt_input['lowest_matched_taxon'].fillna(
                dt_input[col], inplace=True)

        annot_attr = ['rank_spec', 'msms_score', 'libname', 'structure_inchikey', 'structure_inchi', 'structure_smiles', 'structure_molecular_formula', 'adduct',
                      'structure_exact_mass', 'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass',
                      'structure_taxonomy_npclassifier_03class',
                      'query_otol_species', 'lowest_matched_taxon', 'score_taxo', 'score_max_consistency', 'final_score', 'rank_final']

    else:
        annot_attr = ['rank_spec', 'msms_score', 'libname', 'structure_inchikey', 'structure_inchi',
                      'structure_smiles', 'structure_molecular_formula', 'adduct',
                      'structure_exact_mass', 'short_inchikey', 'structure_taxonomy_npclassifier_01pathway',
                      'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class',
                      'organism_name', 'organism_taxonomy_ottid',
                      'organism_taxonomy_01domain', 'organism_taxonomy_02kingdom', 'organism_taxonomy_03phylum',
                      'organism_taxonomy_04class', 'organism_taxonomy_05order', 'organism_taxonomy_06family', 'organism_taxonomy_07tribe',
                      'organism_taxonomy_08genus', 'organism_taxonomy_09species', 'organism_taxonomy_10varietas',
                      'matched_domain', 'matched_kingdom', 'matched_phylum', 'matched_class', 'matched_order',
                      'matched_family', 'matched_tribe', 'matched_genus', 'matched_species', 'score_taxo', 'score_max_consistency', 'final_score', 'rank_final']

    comp_attr = ['component_id', 'structure_taxonomy_npclassifier_01pathway_consensus', 'freq_structure_taxonomy_npclassifier_01pathway',
                 'structure_taxonomy_npclassifier_02superclass_consensus',
                 'freq_structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class_consensus', 'freq_structure_taxonomy_npclassifier_03class']

    col_to_keep = ['feature_id'] + comp_attr + annot_attr

    # We add the min chemo score at this step

    dt_input = dt_input[
        ((dt_input['score_taxo'] >= min_score_taxo_ms1) & (dt_input['score_max_consistency'] >= min_score_chemo_ms1)) | (
            dt_input['libname'] == 'ISDB')]

    dt_output_flat = dt_input[col_to_keep]

    # Cytoscape formatting 

    all_columns = list(dt_input) # Creates list of all column headers
    
    dt_input[all_columns] = dt_input[all_columns].astype(str)

    gb_spec = {c: '|'.join for c in annot_attr}

    for c in comp_attr:
        gb_spec[c] = 'first'

    dt_output_cyto = dt_input.groupby('feature_id').agg(gb_spec)
    dt_output_cyto.reset_index(inplace=True)

    return dt_output_flat, dt_output_cyto


def chembl_id_fetcher(df_input):

    """ This function will get the CHEMBL ids from the structure_inchikey field

    Args:
        dt_isdb_results (dataframe) : a annotation table
        top_to_output (integer): Top N of candidate to keep
    Returns:
        dt_isdb_results_chem_rew (dataframe): a dataframe with the top N annotation ordered by final rank
    """

    # fetching CHEMBL infos

    molecule = new_client.molecule

    inchi_keys = df_input['structure_inchikey'].unique()
    chunks_query = [inchi_keys[x:x+35] for x in range(0, len(inchi_keys), 35)]

    results = []
    bad_keys = []


    for chunk in tqdm(chunks_query):
        try:
            res = molecule.get(list(chunk))
            results.append(res)
        except:
            # Inchi key was not found in ChEMBL
            bad_keys.append(chunk)

    flat_list = [item for sublist in results for item in sublist]

    chembl_df = json_normalize(flat_list)

    # A security in case no chembl results are outputted
    
    if len(chembl_df) > 0: 

        chembl_df = chembl_df[['molecule_chembl_id', 'molecule_structures.standard_inchi_key']]


        df_output = pd.merge(left=df_input, right=chembl_df, left_on='structure_inchikey', right_on='molecule_structures.standard_inchi_key', how = 'left')

        df_output.rename(columns={'molecule_chembl_id': 'structure_chembl_id'}, inplace=True)
        df_output.drop(['molecule_structures.standard_inchi_key'], axis=1, inplace=True)

        return df_output
    else:
        return df_input



def pivot_tabler(df_input, lib_to_keep, minimal_taxo_score, minimal_chemo_score, minimal_total_score, isdb_results_repond_flat_sel_path, pivot_table_results_path):

    if len(lib_to_keep) != 0:
            df_input = df_input[df_input['libname'].str.contains(lib_to_keep)]

    df_input[df_input['score_taxo'] >= minimal_taxo_score]
    df_input[df_input['score_max_consistency'] >= minimal_chemo_score]
    df_input[df_input['final_score'] >= minimal_total_score]


    df_input_sel = df_input[['feature_id', 'component_id', 'structure_taxonomy_npclassifier_01pathway_consensus','structure_taxonomy_npclassifier_02superclass_consensus',
    'structure_taxonomy_npclassifier_03class_consensus', 'msms_score', 'libname',
    'structure_inchikey', 'structure_inchi', 'structure_smiles', 'structure_molecular_formula',
    'adduct', 'structure_exact_mass', 'short_inchikey',
    'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass',
    'structure_taxonomy_npclassifier_03class', 'organism_name', 'organism_taxonomy_ottid',
    'organism_taxonomy_01domain', 'organism_taxonomy_02kingdom', 'organism_taxonomy_03phylum',
    'organism_taxonomy_04class', 'organism_taxonomy_05order', 'organism_taxonomy_06family',
    'organism_taxonomy_07tribe', 'organism_taxonomy_08genus', 'organism_taxonomy_09species',
    'score_taxo', 'score_max_consistency', 'final_score']]

    df_input_sel.to_csv(isdb_results_repond_flat_sel_path, sep='\t', index=None)

    pivot_ui(df_input_sel, outfile_path=pivot_table_results_path)

