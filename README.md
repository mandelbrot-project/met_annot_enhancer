# met_annot_enhancer

A set of script to proceed to metabolite annotation results enhancement using taxonomically and structurally informed processes.

# Installation
## 1. Clone this repo

`git clone https://github.com/mandelbrot-project/met_annot_enhancer.git`


## 2. Set the correct environment


The requested package can be installed by creating a conda environment and activating it.

Use the following line for environment creation 

`conda env create -f environment.yml`

And this one to activate it 

`conda activate met_annot_enhancer_env`

# Requirements 
## A FBMN job on GNPS

For now the scripts work by taking output of Feature Based Molecular Networking jobs.
You can read more about FBMN here (https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/)

## A metadata file containing taxonomical information 

A specific requirement to proceed to taxonomical reweighting of your annotation results is that a metadata file with the taxonomical information of your samples is provided. See below for details.

2 columns of metadata are required for taxonomical reweighting:
  - a *filename* column (with 1 row / MS-file). Example: 2013901_AG_07098.mzML.
        ! The header for this column has to be *filename*.
  - an *organism* column (organism the sample was collected from). Example: *Ambystoma mexicanum*

This file is ideally uploaded when you relize your FBMN job. See here for details on the metadata file format (https://ccms-ucsd.github.io/GNPSDocumentation/metadata/)
 
If you didn't add it you can do this afterward. For this, the metadata file, in the form of a *.tsv* file, should be placed in data_in/yourjobid/metadata_table/*metadata_file*. The name of the file doesn't matter but it should be the only file in the folder.
  
# Usage

## 1. Edit the .yaml file containing the job parameters

For this you can first copy the default.yaml file from [here](https://github.com/mandelbrot-project/met_annot_enhancer/tree/main/configs/default) to 
[there](https://github.com/mandelbrot-project/met_annot_enhancer/tree/main/configs/user_defined). Don't rename the file.

Now you can safely edit the file in configs/user_defined/default.yaml according to your needs. See below for a brief descrpition of each parameters
### Parameters description: 

options:

  - **download_gnps_job: True**
    - set to False it you already downloaded a GNPS FBMN
  - **do_spectral_match: True**
    - will perform MS2 matching using db_file_path .mgf
  - **keep_lowest_taxon: False**
    - for clarity un outputs, just keep the lowest taxon matched
  - **output_plots: True**
    - False if keep_lowest_taxon = True (to change)
 
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
  - **var_one_header: 'SAMPLE_info'**
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
  
## 2. Launch the job

From the home folder of this repository. In the activated conda environment.

`python src/dev/nb.py`


# References

Description and original implementation of the taxonomically informed metabolite annotation process is available here https://doi.org/10.3389/FPLS.2019.01329, associated data have been deposited in the following OSF repository: <https://osf.io/bvs6x/>.
A snapshot of the code at the time of publication is also available at <https://github.com/oolonek/taxo_scorer>.

An R based implementation of the metabolite annotation enhancing process is available here https://github.com/taxonomicallyinformedannotation/tima-r
