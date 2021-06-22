import zipfile
import glob
import os
import sys
import shlex
import subprocess
import yaml

# Loading the parameters from yaml file
print('Loading parameters...')

with open (r'./configs/gnps_param.yaml') as file:    
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    params_list = yaml.load(file, Loader=yaml.FullLoader)

job_id = params_list['job'][0]['job_id']
gnps_job_path = params_list['job'][1]['gnps_job_path']
project_name = params_list['job'][2]['project_name']

base_filename = 'GNPS_output_' + project_name
filename_suffix = 'zip'

path_to_folder = os.path.expanduser(os.path.join(gnps_job_path, base_filename))
path_to_file = os.path.expanduser(os.path.join(gnps_job_path, base_filename + "." + filename_suffix))
base_filename = 'GNPS_output_' + project_name

# Downloading GNPS files
files = glob.glob(gnps_job_path)
for f in files:
    if f == path_to_file:
        os.remove(f)

print('''
Fetching the GNPS job ...
''')

job_url_zip = "https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task="+job_id+"&view=download_cytoscape_data"
print(job_url_zip)

cmd = 'curl -d "" '+job_url_zip+' -o '+path_to_file+ ' --create-dirs'
subprocess.call(shlex.split(cmd))

with zipfile.ZipFile(path_to_file, 'r') as zip_ref:
    zip_ref.extractall(path_to_folder)

# We finally remove the zip file
os.remove(path_to_file)

params_suffix = '.yaml'
with open(os.path.join(path_to_folder, job_id + params_suffix), 'w') as file:  
    documents = yaml.dump(params_list, file)
