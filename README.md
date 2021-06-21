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


## Proceeding to spectral matching and followed by taxonomically and struturally informed scoring 

### 1.  edit the parameters in default.yaml

(detail later on here and in the yaml wich are the important params)

### 2.  launch the job

`python src/dev/met_annot_enhancer_with_params.py`


# R implementation

An R based implementation of the metabolite annotation enhancing process is available here https://gitlab.com/tima5/taxoscorer

# References

Description and original implementation of the taxonomically informed metabolite annotation process is available here https://doi.org/10.3389/FPLS.2019.01329, associated data have been deposited in the following OSF repository: <https://osf.io/bvs6x/>.
A snapshot of the code at the time of publication is also available at <https://github.com/oolonek/taxo_scorer>.
