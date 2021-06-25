# met_annot_enhancer
A set of script to proceed to metabolite annotation results enhancement using taxonomically and structurally informed processes

# Usage

## Requirements 


### 1.  clone this repo

`git clone https://github.com/mandelbrot-project/met_annot_enhancer.git`


### 2.  set the correct environmenet


The requested package can be installed by creating the conda environment and activating it.

Use the following line for environment creation 

`conda env create -f environment.yml`

And this one to activate it 

`conda activate met_annot_enhancer_env`

If you need to update the environment run 

`conda env update --file environment.yml`

## STEP A) Download FBMN from GNPS 

You need first to perform a FBMN in orderto get the component index (ie cluster) of the features, that will be used for the reweighting.
You can download a pre-computed GNPS job using the following commmand line; you just need your job id. Results will be stored in the data_in/yourjobid folder
DO NOT CHANGE THE NAME OF THE FOLDER

`python src/dev/download_from_gnps.py --job_id yourjobidgoeshere`

for help, use:

`python src/dev/download_from_gnps.py -h`

## STEP B) Metadata file required 

2 columns of metadata are required for taxonomical reweighting:
  - a *filename* column (with 1 row / MS-file). Example: 2013901_AG_07098.mzML.
        ! The header for this column has to be *filename*.
  - an *organism* column (organism the sample was collected from). Example: *Ambystoma mexicanum*
 
Add metadata to your GNPS results if you didn't upload it on GNPS when you ran the job, or if the one you uploaded dosens't contain the information needed.
The metadata file, in the form of a .csv, should be placed in data_in/yourjobid/metadata_table/*metadata_file*. The name of the file doen't matter but it should be the only file in the folder.
  
## STEP C) Proceed to spectral matching followed by taxonomically and structurally informed scoring 

### 1.  Copy the default.yaml and edit yccording to your needs

Copy the default.yaml file from [here](https://github.com/mandelbrot-project/met_annot_enhancer/tree/main/configs/default) to 
[there](https://github.com/mandelbrot-project/met_annot_enhancer/tree/main/configs/user_defined). Don't rename the file.

### 2.  edit the parameters in configs/user_defined/default.yaml

Parameters are: 

options:
  - **download_gnps_job: False**
    - set to False it you already downloaded a GNPS FBMN in step A
  - **do_spectral_match: True**
    - will perform MS2 matching using db_file_path .mgf
  - **keep_lowest_taxon: False**
    - for clarity un outputs, just keep the lowest taxon matched
  - **output_plots: True**
    - False
 
paths:
  - **gnps_job_id: 250536f4cb3e4f159e5ef67a3d024fac**
    - The GNPS job id you want to annotate
  - **project_name: your name**
    - The name you want to give to your project, output resulst in data_out/project_name
  - **metadata_path: db_metadata/210523_lotus_dnp_metadata.csv**
    - Path to your spectral library file
  - **db_file_path: db_spectra/LOTUS_DNP_ISDB.mgf**
    - Path to the metadata of the spectral file
  - **adducts_pos_path: data_loc/db_prepared_pos.tsv.gz**
    - Path to the adducts file in pos mode
  - **adducts_neg_path: data_loc/db_prepared_neg.tsv.gz**
    - Path to the adducts file in neg mode

spectral_match_params:
  - **parent_mz_tol: 0.01**
    - the parent mass tolerance to use for spectral matching (in Da)
  - **msms_mz_tol: 0.01**
    - the msms mass tolerance to use for spectral matching (in Da)
  - **min_cos: 0.3**
    - the minimal cosine to use for spectral matching
  - **min_peaks: 8**
    - the minimal matching peaks number to use for spectral matching

repond_params:
  - **Top_N_Sample: 0**
    - Max number of contributors to take into account for taxo reponderation, set to 0 for all biosources where the feature is detected
  - **top_to_output: 1**
    - Top X for final ouput
  - **ppm_tol: 2**
    - ppm tolerance to be used for ms1 match
  - **polarity: 'pos'**
    - ion mode you are working with (pos or neg)
  - **organism_header: 'ATTRIBUTE_Species'**
    - Mandatory: header of your samples' organism in metadata file
  - **sampletype_header: 'SAMPLE_info'**
    - Optional parameter
  - **use_post_taxo: True**
    - Set True if you want to use rank after taxonomical reweighting for consensus chemical class determination, else set to False
  - **top_N_chemical_consistency: 15**
    - if use_post_taxo = True: Top N to use for chemical consistency (annotation not in top N will be discared for component consensus determination)
  - **file_extension: '.mzML'**
    - MS filename extension (or any common pattern in all your MS filenames)
  - **msfile_suffix: ' Peak area'**
    - Suffix to remove in you *quant.csv* table to match your MS filenames in metadata table
  - **min_score_ms1: 5**
    - Minimum taxonomical score (7 = species, 6 = genus, 5 = family, ...)
  
### 2.  launch the job

`python src/dev/met_annot_enhancer.py`


# R implementation

An R based implementation of the metabolite annotation enhancing process is available here https://gitlab.com/tima5/tima-weighter

# References

Description and original implementation of the taxonomically informed metabolite annotation process is available here https://doi.org/10.3389/FPLS.2019.01329, associated data have been deposited in the following OSF repository: <https://osf.io/bvs6x/>.
A snapshot of the code at the time of publication is also available at <https://github.com/oolonek/taxo_scorer>.
