# met_annot_enhancer
A set of script to proceed to metabolite annotation results enhancement using taxonomically and structurally informed processes

# Usage

## Requirements 

The requested package can be installed by creating the conda environment and activating it.

Use the following line for environment creation 

`conda env create -f environment.yml`

And this one to activate it 

`conda activate met_annot_enhancer_env`

If you need to update the environment run 

`conda env update --file environment.yml`


## Launching the script

```
python met_annot_enhancer.py \
gnps_job_id \
gnps_path \
spectral_lib_matcher_result_path \
metadata_path \
met_annot_enhancer_results_path \
top N of desired outputs per feature \
tolerance for MS1 match (in ppm) \
polarity mode (Pos or Neg)
```

Example : 

```
python met_annot_enhancer.py \
56d01c6ccfe143eca5252017202c8fef \
/Users/pma/tmp/Fred_Legendre/ \
/Users/pma/tmp/Fred_Legendre/GNPS_output/spectral_matcher_results_DNP_ISDB.tsv \
/Users/pma/Documents/190602_DNP_TAXcof_CF.tsv \
/Users/pma/tmp/Fred_Legendre/GNPS_output/spectral_matcher_results_DNP_ISDB_repond.tsv \
3 \
5 \
Pos
````


# R implementation

An R based implementation of the metabolite annotation enhancing process is available here https://gitlab.com/tima5/taxoscorer

# References

Description and original implementation of the taxonomically informed metabolite annotation process is available here https://doi.org/10.3389/FPLS.2019.01329, associated data have been deposited in the following OSF repository: <https://osf.io/bvs6x/>.
A snapshot of the code at the time of publication is also available at <https://github.com/oolonek/taxo_scorer>.
