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