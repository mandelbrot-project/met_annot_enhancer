# required libraries

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
import spectral_lib_matcher

def yaml_params_loader(params_list):

    """ Gets individuals params from a yaml file
    Args:
            a params_list(list): a Glist of parameters loaded from a list
    Returns:
            individual parameters
    """

    download_gnps_job = params_list['options']['download_gnps_job']
    do_spectral_match = params_list['options']['do_spectral_match']
    do_taxo_resolving = params_list['options']['do_taxo_resolving']
    keep_lowest_taxon = params_list['options']['keep_lowest_taxon']
    output_plots = params_list['options']['output_plots']


    gnps_job_id = params_list['paths']['gnps_job_id']
    input_folder = params_list['paths']['input_folder']
    project_name = params_list['paths']['project_name']
    output_folder = params_list['paths']['output_folder']
    metadata_path = params_list['paths']['metadata_path']
    db_file_path = params_list['paths']['db_file_path']
    adducts_pos_path = params_list['paths']['adducts_pos_path']
    adducts_neg_path = params_list['paths']['adducts_neg_path']

    parent_mz_tol = params_list['spectral_match_params']['parent_mz_tol']
    msms_mz_tol = params_list['spectral_match_params']['msms_mz_tol']
    min_cos = params_list['spectral_match_params']['min_cos']
    min_peaks = params_list['spectral_match_params']['min_peaks']

    Top_N_Sample = params_list['repond_params']['Top_N_Sample']
    top_to_output= params_list['repond_params']['top_to_output']
    ppm_tol = params_list['repond_params']['ppm_tol']
    polarity = params_list['repond_params']['polarity']
    organism_header = params_list['repond_params']['organism_header']
    sampletype_header = params_list['repond_params']['sampletype_header']
    use_post_taxo = params_list['repond_params']['use_post_taxo'] # Note we dont use this variable see why !
    top_N_chemical_consistency = params_list['repond_params']['top_N_chemical_consistency']
    file_extension = params_list['repond_params']['file_extension']
    msfile_suffix = params_list['repond_params']['msfile_suffix']
    min_score_taxo_ms1 = params_list['repond_params']['min_score_taxo_ms1']
    min_score_chemo_ms1 = params_list['repond_params']['min_score_chemo_ms1']

    drop_blanks = params_list['plotting_params']['drop_blanks']
    multi_plot = params_list['plotting_params']['multi_plot']

    return gnps_job_id,input_folder,project_name = params_list['paths'][2]['project_name']
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


def gnps_job_fetcher(gnps_job_id, input_folder):

    """Fetch a GNPS job and saves its to a given folder
    Args:
            gnps_job_id (str): a GNPS Job id
            input_folder (str): the folder where the job should be kept
    Returns:
        nothing
    """
    path_to_folder = os.path.expanduser(os.path.join(input_folder , gnps_job_id))
    path_to_file = os.path.expanduser(os.path.join(input_folder , gnps_job_id + '.zip'))

    print('''
    Fetching the GNPS job: '''
    + gnps_job_id
    )
    # or &view=download_clustered_spectra or download_cytoscape_data (check when to use which and optionalize)
    job_url_zip = "https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task="+gnps_job_id+"&view=download_cytoscape_data"
    
    cmd = 'curl -d "" '+job_url_zip+' -o '+path_to_file+ ' --create-dirs'
    subprocess.call(shlex.split(cmd))

    with zipfile.ZipFile(path_to_file, 'r') as zip_ref:
        zip_ref.extractall(path_to_folder)

    # We finally remove the zip file
    os.remove(path_to_file)

    print('''
    Job successfully downloaded: results are in: '''
    + path_to_folder
    )