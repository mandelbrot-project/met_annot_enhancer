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

`python src/dev/download_from_gnps.py --job_id yourjobidgoeshere`

for help, use:

`python src/dev/download_from_gnps.py -h`

## STEP B) Proceeding to spectral matching and followed by taxonomically and struturally informed scoring 

### 1.  edit the parameters in default.yaml

This file is located in the configs folder [here](https://github.com/mandelbrot-project/met_annot_enhancer/blob/194dcde9383f63549241694f2b2ac85635a6f15f/configs/default.yaml)

(detail later on here and in the yaml wich are the important params)

### 2.  launch the job

`python src/dev/met_annot_enhancer_with_params.py`


# R implementation

An R based implementation of the metabolite annotation enhancing process is available here https://gitlab.com/tima5/tima-weighter

# References

Description and original implementation of the taxonomically informed metabolite annotation process is available here https://doi.org/10.3389/FPLS.2019.01329, associated data have been deposited in the following OSF repository: <https://osf.io/bvs6x/>.
A snapshot of the code at the time of publication is also available at <https://github.com/oolonek/taxo_scorer>.
