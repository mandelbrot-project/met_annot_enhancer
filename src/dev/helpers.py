# required libraries

import pandas as pd
import numpy as np
import zipfile
import glob
import os
import sys
import shlex
import subprocess
from opentree import OT
import json
from pandas import json_normalize


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


    paths_dic['quantification_table_reformatted_path'] = os.path.join(path_to_gnps_folder,'quantification_table_reformatted','')

    output_folder = params_list['paths']['output_folder']
    project_name = params_list['paths']['project_name']

    paths_dic['path_to_results_folders'] =  os.path.expanduser(os.path.join(output_folder, project_name +'/'))
    if not os.path.exists(paths_dic['path_to_results_folders']):
        os.makedirs(paths_dic['path_to_results_folders'])
    
    paths_dic['clusterinfo_summary_path'] = os.path.join(path_to_gnps_folder,'clusterinfo_summary','')


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