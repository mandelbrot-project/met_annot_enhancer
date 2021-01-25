# %% Load packages

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





# %% defininbg display options

pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 100)

# We deactivate the iloc warning see https://stackoverflow.com/a/20627316
pd.options.mode.chained_assignment = None  # default='warn'


# %% Choose line by line or classic (1 species) reweighting

## //TODO Maybe remove this option and systematically apply it ?
 
Run_line_x_line = True 
Top_N_Sample = 5

# species_bio = 'Arabidopsis thaliana'
# genus_bio = 'Arabidopsis'
# family_bio = 'Brassicaceae'
# order_bio = 'Brassicales'
# class_bio = 'Magnoliopsida'
# phylum_bio = 'Tracheophyta'
# kingdom_bio = 'Plantae'


# %% Defining the paths

# # GNPS FBMN Job ID
# job_id = "56d01c6ccfe143eca5252017202c8fef"
# # toy set : f840934ade8744f9a6b4438804af8287
# # LP: 75c64084d1e749748d5a29f671944231

# # this is the path for the ISDB resukts file (topX)
# #isdb_results_path = "J:/COMMON FASIE-FATHO/PF_project/Toy_Set/Toy_set_top_50_pos.out"
# isdb_results_path = "/Users/pma/tmp/Fred_Legendre/GNPS_output/spectral_matcher_results_DNP_ISDB.tsv"
# # isdb_results_path = "/Users/pma/tmp/bafu_ecometabo/FBMN_bafu_ecometabo_pos/FBMN_bafu_ecometabo_pos_msmatched_ISDB_DNP.out"


# # this is the path for the DNP or OPENNPDB datatable file
# metadata_path = "/Users/pma/Documents/190602_DNP_TAXcof_CF.tsv"

# #Path to weighed annotation result of ISDB
# output_weighed_ISDB_path = "/Users/pma/tmp/Fred_Legendre/GNPS_output/spectral_matcher_results_DNP_ISDB_repond.tsv"


# # Path for the GNPS job export dir 

# gnps_job_path = "/Users/pma/tmp/Fred_Legendre/"


# # Set True if you want to use rank after taxonomical reweighting for consensus chemical class determination
use_post_taxo = True

# # MS filename extension (a common pattern in all your filenames)
file_extension = '.mzML'

# # Set True if you want to use rank after taxonomical reweighting for consensus chemical class determination
top_N_chemical_consistency = 30

# top_to_output = 3

# ppm_tol = 5

# #polarity = 'Neg'
# polarity = 'Pos'


# defining the command line arguments
try:
    job_id = sys.argv[1]
    gnps_job_path = sys.argv[2]
    isdb_results_path = sys.argv[3]
    metadata_path = sys.argv[4]
    output_weighed_ISDB_path = sys.argv[5]
    top_to_output = sys.argv[6]
    ppm_tol = sys.argv[7]
    polarity = sys.argv[8]

    print('Proceeding to metabolite annotation enhancement on GNPS job id: '
          + job_id + '.\n'
          + 'The GNPS files are downloaded locally in the following folder:' 
          + gnps_job_path + '.\n'
          + 'The metabolite annotation enhancement is done using the following parameters: \n' 
          + '   - spectral match results: ' + str(isdb_results_path) + '\n'
          + '   - metadata file: ' + str(metadata_path) + '\n'
          + '   - top hits to output: ' + str(top_to_output) + '\n'
          + '   - parent mass tolerance: ' + str(ppm_tol) + '\n'
          + '   - polarity mode: ' + str(polarity) + '\n'
          + 'Results will be outputed in the follwoing file: ' + output_weighed_ISDB_path)
except:
    print(
        '''Please add input and output file path as first and second argument, InChI column header as third argument and finally the number of cpus you want to use.
        Example :
        python spectral_lib_matcher.py /Users/pma/tmp/Lena_metabo_local/FBMN_metabo_lena/spectra/fbmn_lena_metabo_specs_ms.mgf /Users/pma/tmp/New_DNP_full_pos.mgf 0.01 0.01 0.2 6 /Users/pma/tmp/lena_matched.out''')

# python met_annot_enhancer.py 
job_id = '7f1259a161974b9fa4215b1f2a6dca5e'
gnps_job_path = '/Users/pma/tmp/bafu_ecometabo/'
isdb_results_path = '/Users/pma/tmp/bafu_ecometabo/GNPS_output/bafu_ecometabo_spectral_match_results.tsv'
metadata_path = '/Users/pma/Documents/190602_DNP_TAXcof_CF.tsv'
output_weighed_ISDB_path = '/Users/pma/tmp/bafu_ecometabo/GNPS_output/bafu_ecometabo_spectral_match_results_repond.tsv'
top_to_output = 3
ppm_tol = 5
polarity = 'Pos'
organism_header = 'sample_type'


# python met_annot_enhancer.py \
# 56d01c6ccfe143eca5252017202c8fef \
# /Users/pma/tmp/Fred_Legendre/ \
# /Users/pma/tmp/Fred_Legendre/GNPS_output/spectral_matcher_results_DNP_ISDB.tsv \
# /Users/pma/Documents/190602_DNP_TAXcof_CF.tsv \
# /Users/pma/tmp/Fred_Legendre/GNPS_output/spectral_matcher_results_DNP_ISDB_repond.tsv \
# 3 \
# 5 \
# Pos

# timer is started
start_time = time.time()


# %% Import MS1 list

if polarity == 'Pos':
    adducts_df = pd.read_csv(
        '../data_loc/db_pos.tsv.gz', compression='gzip', sep='\t')
else:
    adducts_df = pd.read_csv(
        '../data_loc/db_neg.tsv.gz', compression='gzip', sep='\t')

adducts_df['min'] = adducts_df['adduct_mass'] - \
    int(ppm_tol) * (adducts_df['adduct_mass'] / 1000000)
adducts_df['max'] = adducts_df['adduct_mass'] + \
    int(ppm_tol) * (adducts_df['adduct_mass'] / 1000000)


# %% Downloading GNPS files

# files = glob.glob(gnps_job_path)
# for f in files:
#    os.remove(f)

base_filename = 'GNPS_output'
filename_suffix = 'zip'

path_to_folder = os.path.join(gnps_job_path, base_filename)
path_to_file = os.path.join(gnps_job_path, base_filename + "." + filename_suffix)


job_url_zip = "https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task="+job_id+"&view=download_cytoscape_data"

cmd = 'curl -d "" '+job_url_zip+' -o '+path_to_file
subprocess.call(shlex.split(cmd))

with zipfile.ZipFile(path_to_file, 'r') as zip_ref:
    zip_ref.extractall(path_to_folder)

# We finally remove the zip file
cmd = 'rm '+ path_to_file
subprocess.call(shlex.split(cmd))

# %% Spectral matching stage

# Check wether we can import the spectral matching stage here directly

# %% Loading the files

# dt_isdb_results = pd.read_csv(isdb_results_path,
#                               sep='\t',
#                               usecols=['cluster index', 'componentindex', 'parent mass', 'Spectral_Score_DNP',
#                                        'IS_libname', 'IS_match_mzerror', 'InChIKey_DNP'],
#                               error_bad_lines=False, low_memory=True)

## matchms alternative 

dt_isdb_results = pd.read_csv(isdb_results_path,
                              sep='\t',
                              usecols=['msms_score', 'feature_id', 'reference_id', 'inchikey'],
                              error_bad_lines=False, low_memory=True)


## we add a fixed libname (to be changed later on) 

dt_isdb_results['libname'] = 'DNP_ISDB'

dt_isdb_results.rename(columns={'componentindex': 'component_id',
                                'parent mass': 'mz', 'msms_score': 'score_input'}, inplace=True)


## In fact we can directly start here
## we get the networks info (cluster id, component index and parent mass form the downloaded dolder)

clusterinfo_summary_path = os.path.join(path_to_folder,'clusterinfo_summary','')

clusterinfo_summary = pd.read_csv(clusterinfo_summary_path + str(os.listdir(clusterinfo_summary_path)[0]),
                                  sep='\t',
                                  usecols=['cluster index', 'componentindex', 'parent mass'],
                                  error_bad_lines=False, low_memory=True)

clusterinfo_summary.rename(columns={'cluster index': 'feature_id', 'componentindex': 'component_id',
                                'parent mass': 'mz', 'msms_score': 'score_input'}, inplace=True)


# ## we now merge this back with the isdb matched results 
dt_isdb_results = pd.merge(dt_isdb_results, clusterinfo_summary, on='feature_id')

# ## we drop the duplicated ['cluster index'] column
#dt_isdb_results.drop(columns=['cluster index'], inplace=True)

# ## we return a short_inchikey column

# dt_isdb_results['short_inchikey'] = dt_isdb_results['inchikey'].str.split("-", n=1, expand=True)[0]


dt_metadata = pd.read_csv(metadata_path,
                          sep='\t', error_bad_lines=False, low_memory=True)

cluster_count = dt_isdb_results.drop_duplicates(
    subset=['feature_id', 'component_id']).groupby("component_id").count()
cluster_count = cluster_count[['feature_id']].rename(
    columns={'feature_id': 'ci_count'}).reset_index()

print('Number of features: ' + str(len(dt_isdb_results)))
print('Number of annotated features: ' +
      str(dt_isdb_results['inchikey'].count()))


# %%
# Now we directly do the MS1 matching stage on the cluster_summary. No need to have MS2 annotations

super_df = []

for i, row in tqdm(clusterinfo_summary.iterrows(), total=clusterinfo_summary.shape[0]):

    par_mass = clusterinfo_summary.loc[i, 'mz']

    df_0 = clusterinfo_summary.loc[[i], ['feature_id', 'mz']]

    df_1 = adducts_df[(adducts_df['min'] <= par_mass) & (adducts_df['max'] >= par_mass)]
    
    df_1['key'] = i
    df_1.drop(['min', 'max'], axis=1, inplace=True)

    df_tot = pd.merge(df_0, df_1, left_index=True, right_on='key', how='left')
    super_df.append(df_tot)

df_MS1 = pd.concat(super_df, axis=0)
del super_df, adducts_df

df_MS1 = df_MS1.drop(['key'], axis=1).drop_duplicates(
    subset=['feature_id', 'adduct'])

print('MS1 annotation done')


# %%
df_meta_2 = dt_metadata[['IK_DNP', 'Accurate_Mass_DNP']]
df_meta_2.rename(columns={'IK_DNP': 'inchikey'}, inplace=True)

df_meta_2 = df_meta_2.dropna(subset=['Accurate_Mass_DNP'])
df_meta_2 = df_meta_2.drop_duplicates(
    subset=['inchikey', 'Accurate_Mass_DNP'])

df_meta_2 = df_meta_2.round({'Accurate_Mass_DNP': 5})
df_MS1 = df_MS1.round({'exact_mass': 5})

df_MS1_merge = pd.merge(df_MS1, df_meta_2, left_on='exact_mass',
                        right_on='Accurate_Mass_DNP', how='left')
df_MS1_merge = df_MS1_merge.dropna(subset=['inchikey'])

df_MS1_merge['match_mzerror_MS1'] = df_MS1_merge['mz'] - \
    df_MS1_merge['adduct_mass']
df_MS1_merge = df_MS1_merge.round({'match_mzerror_MS1': 5}).astype({
    'match_mzerror_MS1': 'str'})

df_MS1_merge = df_MS1_merge.drop(
    ['Accurate_Mass_DNP', 'adduct_mass', 'mz', 'exact_mass'], axis=1)
df_MS1_merge['score_input'] = 0

#df_MS1_merge = df_MS1_merge.astype({'score_input': 'str'})

# the line below trhrow an error on my side 
# TypeError: groupby() got an unexpected keyword argument 'dropna'
# It should normally be fixed in pandas 1.2.0 As Arnaud which version he is using (probably an older one)
# see https://github.com/pandas-dev/pandas/issues/37323 for a quick and dirty fix
# Actually I even think no dropna is needed
df_MS1_merge_gb = df_MS1_merge.groupby(['feature_id'], dropna=False).agg('|'.join)

df_MS1_merge_gb = df_MS1_merge.groupby(['feature_id']).agg('|'.join)
df_MS1_merge_gb.reset_index(inplace=True)

# %%

dt_isdb_results = pd.concat([dt_isdb_results, df_MS1_merge])
dt_isdb_results.info()




# dt_isdb_results = pd.merge(
#     dt_isdb_results, df_MS1_merge_gb, on='feature_id', how='left')

# dt_isdb_results["inchikey"] = dt_isdb_results["inchikey"].str.cat(
#     dt_isdb_results["inchikey_MS1"], sep="|")
# dt_isdb_results["score_input"] = dt_isdb_results["score_input"].str.cat(
#     dt_isdb_results["score_input_MS1"], sep="|")
# dt_isdb_results["match_mzerror"] = dt_isdb_results["match_mzerror"].str.cat(
#     dt_isdb_results["match_mzerror_MS1"], sep="|")
# dt_isdb_results["libname"] = dt_isdb_results["libname"].str.cat(
#     dt_isdb_results["adduct"], sep="|")

# dt_isdb_results['inchikey'].fillna(
#     dt_isdb_results['inchikey_MS1'], inplace=True)
# dt_isdb_results['score_input'].fillna(
#     dt_isdb_results['score_input_MS1'], inplace=True)
# dt_isdb_results['match_mzerror'].fillna(
#     dt_isdb_results['match_mzerror_MS1'], inplace=True)
# dt_isdb_results['libname'].fillna(dt_isdb_results['adduct'], inplace=True)

# dt_isdb_results = dt_isdb_results.drop(
#     ['inchikey_MS1', 'score_input_MS1', 'mz', 'match_mzerror_MS1', 'adduct'], axis=1)


# print('Number of annotated features after MS1: ' +
#       str(dt_isdb_results['inchikey'].count()))


print('Number of annotated features after MS1: ' +
      str(len(df_MS1_merge['feature_id'].unique())))


len(clusterinfo_summary['feature_id'].unique())

# %% Here we want to split the unknown numbers of results (in the score and IK columns) into multiple columns
# These stages are not mandatory anymore since we have a long file as output of matchms

# dt_isdb_results = dt_isdb_results.join(
#     dt_isdb_results['score_input'].str.split(
#         '|', expand=True).add_prefix('score_input_')
# ).join(dt_isdb_results['inchikey'].str.split('|', expand=True).add_prefix('inchikey_')
#        ).join(dt_isdb_results['match_mzerror'].str.split('|', expand=True).add_prefix('match_mzerror_')
#               ).join(dt_isdb_results['libname'].str.split('|', expand=True).add_prefix('libname_')
#                      )

# dt_isdb_results.drop(
#     columns=['score_input', 'inchikey', 'libname', 'match_mzerror'], inplace=True)


# %%
# Now we melt the whole stuff. For this we use the wide to long function https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.wide_to_long.html

# dt_isdb_results = pd.wide_to_long(dt_isdb_results, stubnames=[
#                                   "score_input", "inchikey", 'match_mzerror', 'libname'], i='feature_id', j="rank_spec", sep='_')
# dt_isdb_results.reset_index(inplace=True)
# dt_isdb_results.dropna(subset=['inchikey'], axis=0, inplace=True)

print('Total number of MSMS annotations: ' + str(len(dt_isdb_results)))

# %%
# Rank annotations based on the spectral score

dt_isdb_results["score_input"] = pd.to_numeric(
    dt_isdb_results["score_input"], downcast="float")
dt_isdb_results['rank_spec'] = dt_isdb_results[['feature_id', 'score_input']].groupby(
    'feature_id')['score_input'].rank(method='dense', ascending=False)

dt_isdb_results.info()

# %%
# Joining the DNP metadata

# we start by outputing the SIK for the ISDB output

dt_isdb_results['short_inchikey'] = dt_isdb_results.inchikey.str.split(
    "-", expand=True)[0]
dt_isdb_results.reset_index(inplace=True)

# now we merge with the DNP metadata after selection of our columns of interest

cols_to_use = ['CRC_Number_DNP', 'IK_DNP', 'InChI_DNP',
               'Molecule_Name_DNP', 'Molecule_Formula_DNP',
               'Accurate_Mass_DNP', 'Short_IK_DNP', 'Kingdom_cof_DNP', 'Phylum_cof_DNP',
               'Class_cof_DNP', 'Order_cof_DNP', 'Family_cof_DNP', 'Genus_cof_DNP',
               'Species_cof_DNP', 'ClassyFy_Status_DNP',
               'Kingdom_cf_DNP', 'Superclass_cf_DNP', 'Class_cf_DNP',
               'Subclass_cf_DNP', 'Parent_Level_1_cf_DNP']

dt_isdb_results.dropna(subset=['short_inchikey'], inplace=True)
dt_isdb_results = pd.merge(
    left=dt_isdb_results, right=dt_metadata[cols_to_use], left_on='short_inchikey', right_on='Short_IK_DNP', how='outer')
dt_isdb_results.dropna(subset=['feature_id'], inplace=True)


# dt_isdb_results = dt_isdb_results.astype({'feature_id' : 'int64'})


print('Total number of annotations with unique Biosource/line: ' +
      str(len(dt_isdb_results)))

# %% Resolving the taxon information from the GNPS metadata file

# the metadata table path is generated from the base bath to the GNPS results folder
metadata_table_path = os.path.join(path_to_folder,'metadata_table','')

# the metadata table is loaded using the organism column specified before

samples_metadata = pd.read_csv(metadata_table_path + str(os.listdir(metadata_table_path)[0]), sep='\t',
                                   usecols=['filename', organism_header])

# Now we want to get the taxonomic information for each of the samples
# so we want to extract the species information from the metadata file
# %%
samples_metadata[organism_header].dropna(inplace = True)
samples_metadata[organism_header] = samples_metadata[organism_header].str.lower()
species = samples_metadata[organism_header].unique()
len_species = len(species)

print("%s unique species have been selected from the metadata table." % len_species )
# %%

species_tnrs_matched = OT.tnrs_match(species, context_name=None, do_approximate_matching=True, include_suppressed=False)



# %%

with open('species.json', 'w') as out:
    sf = json.dumps(species_tnrs_matched.response_dict, indent=2, sort_keys=True)
    out.write('{}\n'.format(sf))

# %%
with open("species.json") as tmpfile:
        jsondic = json.loads(tmpfile.read())

json_normalize(jsondic)
# %%

df_species_tnrs_matched = json_normalize(jsondic,
               record_path=['results', 'matches']
               )



df_species_tnrs_unmatched = json_normalize(jsondic,
               record_path=['unmatched_names']
               )
# %%

df_species_tnrs_matched.info()


# %%
len(df_species_tnrs_matched['taxon.ott_id'].unique())
# %%


# We then want to match with the accepted name instead of the synonym in case both are present. 
# We thus order by matched_name and then by is_synonym status prior to returning the first row.

df_species_tnrs_matched.sort_values(['search_string', 'is_synonym'], axis = 0, inplace = True)
df_species_tnrs_matched_unique = df_species_tnrs_matched.drop_duplicates('search_string', keep = 'first')

# both df are finally merged
merged_df = pd.merge(samples_metadata, df_species_tnrs_matched_unique, how='left', left_on=organism_header, right_on='search_string', indicator=True)


# %%
#Now we want to retrieve the upper taxa lineage for all the samples

# Firsst we retrieve a list of unique ott_ids

# Here when checking the columns datatype we observe that the ott_ids are as float.
# We need to keep them as int
# displaying the datatypes 
display(merged_df.dtypes) 

# converting 'ott_ids' from float to int (check the astype('Int64') whic will work while the astype('int') won't see https://stackoverflow.com/a/54194908)
merged_df['taxon.ott_id'] = merged_df['taxon.ott_id'].astype('Int64')
  

# However, we then need to put them back to 
merged_df['taxon.ott_id']
ott_list = list(merged_df['taxon.ott_id'].dropna().astype('int'))

#ott_list = ott_list[0:10]

# %%

taxon_info = []

for i in ott_list:
    query = OT.taxon_info(i, include_lineage=True)
    taxon_info.append(query)

# %%


tl = []

for i in taxon_info:
    with open('taxon_info.json', 'w') as out:
        tl.append(i.response_dict)
        yo = json.dumps(tl)
        out.write('{}\n'.format(yo))

# %%

with open("taxon_info.json") as tmpfile:
        jsondic = json.loads(tmpfile.read())

df = json_normalize(jsondic)


# %%

df_tax_lineage = json_normalize(jsondic,
               record_path=['lineage'],
               meta = ['ott_id', 'unique_name'],
               record_prefix='sub_',
               errors='ignore'
               )
# %%
# This keeps the last occurence of each ott_id / sub_rank grouping https://stackoverflow.com/a/41886945

df_tax_lineage_filtered = df_tax_lineage.groupby(['ott_id', 'sub_rank'], as_index=False).last()
# %%
#Here we pivot long to wide to get the taxonomy

df_tax_lineage_filtered_flat = df_tax_lineage_filtered.pivot(index='ott_id', columns='sub_rank', values='sub_name')

# %%
# Here we actually also want the lowertaxon (species usually) name

df_tax_lineage_filtered_flat = pd.merge(df_tax_lineage_filtered_flat, df_tax_lineage_filtered[['ott_id', 'unique_name']], how='left', on='ott_id', )

#Despite the left join ott_id are duplicated 

df_tax_lineage_filtered_flat.drop_duplicates(subset = ['ott_id', 'unique_name'], inplace = True)

# %%
# we keep the fields of interest

df_tax_lineage_filtered_flat[['ott_id', 'kingdom', 'phylum',
                              'class', 'order', 'family', 'genus', 'unique_name']]



# We now rename our columns of interest

renaming_dict = {'kingdom': 'query_otol_kingdom',
                 'phylum': 'query_otol_phylum',
                 'class': 'query_otol_class',
                 'order': 'query_otol_order',
                 'family': 'query_otol_family',
                 'genus': 'query_otol_genus',
                 'unique_name': 'query_otol_species'}


df_tax_lineage_filtered_flat.rename(columns=renaming_dict, inplace=True)

# We select columns of interest 

cols_to_keep = ['ott_id',
                'query_otol_kingdom',
                'query_otol_phylum',
                'query_otol_class',
                'query_otol_order',
                'query_otol_family',
                'query_otol_genus',
                'query_otol_species']

df_tax_lineage_filtered_flat = df_tax_lineage_filtered_flat[cols_to_keep]


# We merge this back with the samplemetadata

samples_metadata = pd.merge(merged_df, df_tax_lineage_filtered_flat, how='left', left_on='taxon.ott_id', right_on='ott_id' )




# %% Extracting biosource / feature for line by line


quantification_table_reformatted_path = os.path.join(path_to_folder,'quantification_table_reformatted','')

metadata_table_path = os.path.join(path_to_folder,'metadata_table','')


if Run_line_x_line == True:

    feature_intensity = pd.read_csv(quantification_table_reformatted_path + str(
        os.listdir(quantification_table_reformatted_path)[0]), sep=',')
    feature_intensity.rename(columns={'row ID': 'row_ID'}, inplace=True)
    feature_intensity.set_index('row_ID', inplace=True)
    feature_intensity = feature_intensity.filter(
        regex=file_extension + '|row_ID')
    feature_intensity = feature_intensity.where(feature_intensity.apply(
        lambda x: x.isin(x.nlargest(Top_N_Sample)), axis=1), 0)  # top N here
    feature_intensity.columns = feature_intensity.columns.str.replace(
        r' Peak area', '')
    feature_intensity = feature_intensity.transpose()
    feature_intensity.index.name = 'MS_filename'
    feature_intensity = feature_intensity.transpose()
    # Samples_metadata = pd.read_csv(metadata_table_path + str(os.listdir(metadata_table_path)[0]), sep='\t',
    #                                # usecols=['filename','ATTRIBUTE_phylum_cof', 'ATTRIBUTE_kingdom_cof',  'ATTRIBUTE_class_cof', 'ATTRIBUTE_order_cof', 'ATTRIBUTE_family_cof', 'ATTRIBUTE_genus_cof', 'ATTRIBUTE_species_cof'])
    #                                #    usecols=['filename', 'ATTRIBUTE_Phylum', 'ATTRIBUTE_Kingdom',  'ATTRIBUTE_Class', 'ATTRIBUTE_Order', 'ATTRIBUTE_Family', 'ATTRIBUTE_Genus', 'ATTRIBUTE_Species'])
    #                                usecols=['filename',
    #                                         'query_otol_kingdom',
    #                                         'query_otol_phylum',
    #                                         'query_otol_class',
    #                                         'query_otol_order',
    #                                         'query_otol_family',
    #                                         'query_otol_genus',
    #                                         'query_otol_species'])
    res = feature_intensity[feature_intensity != 0].stack()
    df_res = res.to_frame().reset_index()
    df_merged = pd.merge(df_res, samples_metadata, left_on='MS_filename',
                         right_on='filename', how='left').drop([0, 'MS_filename', 'filename'], axis=1)
    df_merged = df_merged.groupby('row_ID').agg(lambda x: list(x))
    df_merged.reset_index(inplace=True)


# %%

# Here we will add three columns (even for the simple repond this way it will be close to the multiple species repond)
# these line will need to be defined as function arguments

if Run_line_x_line == True:
    dt_isdb_results = pd.merge(
        dt_isdb_results, df_merged, left_on='feature_id', right_on='row_ID', how='left')
else:
       dt_isdb_results['query_otol_species'] = species_bio
       dt_isdb_results['query_otol_genus'] = genus_bio
       dt_isdb_results['query_otol_family'] = family_bio
       dt_isdb_results['query_otol_order'] = order_bio
       dt_isdb_results['query_otol_class'] = class_bio
       dt_isdb_results['query_otol_phylum'] = phylum_bio
       dt_isdb_results['query_otol_kingdom'] = kingdom_bio
       
#%% Taxonomical Reweighting

cols_ref = ['Kingdom_cof_DNP', 'Phylum_cof_DNP',  'Class_cof_DNP', 'Order_cof_DNP', 'Family_cof_DNP', 'Genus_cof_DNP', 'Species_cof_DNP']
#cols_att = ['ATTRIBUTE_Kingdom', 'ATTRIBUTE_Phylum',  'ATTRIBUTE_Class', 'ATTRIBUTE_Order', 'ATTRIBUTE_Family', 'ATTRIBUTE_Genus', 'ATTRIBUTE_Species']
#cols_att = ['ATTRIBUTE_kingdom_cof', 'ATTRIBUTE_phylum_cof', 'ATTRIBUTE_class_cof', 'ATTRIBUTE_order_cof', 'ATTRIBUTE_family_cof', 'ATTRIBUTE_genus_cof', 'ATTRIBUTE_species_cof']
cols_att = ['query_otol_kingdom',
            'query_otol_phylum',
            'query_otol_class',
            'query_otol_order',
            'query_otol_family',
            'query_otol_genus',
            'query_otol_species']
cols_match = ['matched_kingdom', 'matched_phylum', 'matched_class', 'matched_order', 'matched_family', 'matched_genus', 'matched_species']

col_prev = None
if Run_line_x_line == True:
       for col_ref, col_att, col_match in zip(cols_ref, cols_att, cols_match):
              dt_isdb_results[col_ref].fillna('Unknown', inplace=True)
              dt_isdb_results[col_ref] = dt_isdb_results[col_ref].apply(lambda x: [x])
            #   print(dt_isdb_results[col_ref])
            #   print(dt_isdb_results[col_att])
            #   #dt_isdb_results[col_att] = dt_isdb_results[col_att].apply(lambda x: tuple(x.split(","))) #if original biosource as string
            #   dt_isdb_results[col_ref] = dt_isdb_results[col_ref].apply(lambda x: list(x))
            #   dt_isdb_results[col_att] = dt_isdb_results[col_att].apply(lambda x: list(x))
              #dt_isdb_results[col_att] = dt_isdb_results[col_att].tolist()
              dt_isdb_results[col_match] = [list(set(a).intersection(set(b))) for a, b in zip(dt_isdb_results[col_ref], dt_isdb_results[col_att])] # Allows to compare 2 lists
              dt_isdb_results[col_match] = dt_isdb_results[col_match].apply(lambda y: np.nan if len(y)==0 else y)
              if col_prev != None:
                     dt_isdb_results[col_match].where(dt_isdb_results[col_prev].notnull(), np.nan)
              col_prev = col_match

else:
       for col_ref, col_att, col_match in zip(cols_ref, cols_att, cols_match):
              dt_isdb_results[col_ref].fillna('Unknown', inplace=True)
              dt_isdb_results[col_match] = np.where((dt_isdb_results[col_ref] == dt_isdb_results[col_att]), dt_isdb_results[col_att], np.nan)
              if col_prev != None:
                     dt_isdb_results[col_match].where(dt_isdb_results[col_prev].notnull(), np.nan)
              col_prev = col_match


dt_isdb_results['score_taxo'] = dt_isdb_results[cols_match].count(axis=1)


# [list(set(a).intersection(set(b))) for a, b in zip(dt_isdb_results[col_ref], dt_isdb_results[col_att])]

# a = 'Lotus corniculatus'
# b = ['Arrhenatherum elatus', 'Lotus corniculatus', 'Lotus corniculatus']
# list(set(a).intersection(set(b)))


# dt_isdb_results[col_ref].apply(lambda x: [x])

# # We get a AttributeError: 'tuple' object has no attribute 'split'
# #for dt_isdb_results[col_att].apply(lambda x: tuple(x.split(",")))#



# for a, b in zip(dt_isdb_results[col_ref], dt_isdb_results[col_att]) :
#     print(a, b)
# # type(dt_isdb_results[col_att][0])
# # dt_isdb_results[col_att][0][0]

# # boh = zip(dt_isdb_results[col_att][0])

# # unlist(boh)

# # dt_isdb_results[col_att].apply(lambda x: tuple(x.split(",")))
# # dt_isdb_results[col_att].apply(lambda x: dt_isdb_results[col_att][][x])


# dt_isdb_results[cols_att]
# dt_isdb_results[cols_ref]


# # Put your dataframe here
# df = pd.DataFrame({'A':[1,2], 'B':[(1,2), (3,4)]})  

# print("Original Dataset")
# print(df)

# start = time.time()
# df[['B1','B2']] = pd.DataFrame(df['B'].tolist(),index=df.index)
# print("Method 1")
# print("Time elapsed :" + str(time.time()-start))
# print(df)

# pd.DataFrame(df['B'].tolist(),index=df.index)

# pd.DataFrame(dt_isdb_results[col_att].tolist())


# pd.DataFrame(dt_isdb_results['ATTRIBUTE_Species'].tolist().tolist())

# dt_isdb_results['ATTRIBUTE_Species'].apply(lambda x: tuple(x.split(",")))



# %%
# Filter out MS1 annotations without a reweighting at the order level at least

if polarity == 'Pos':
    dt_isdb_results = dt_isdb_results[(dt_isdb_results['score_taxo'] >= 4) | (
        dt_isdb_results['libname'] == 'DNP_ISDB')]
else:
    dt_isdb_results = dt_isdb_results[(dt_isdb_results['score_taxo'] >= 4) | (
        dt_isdb_results['libname'] == 'DNP_ISDB')]

print('Total number of annotations after filtering MS1 annotations not reweighted at order level: ' +
      str(len(dt_isdb_results)))
print('Number of annotations reweighted at the kingdom level: ' +
      str(dt_isdb_results['matched_kingdom'].count()))
print('Number of annotations reweighted at the phylum level: ' +
      str(dt_isdb_results['matched_phylum'].count()))
print('Number of annotations reweighted at the class level: ' +
      str(dt_isdb_results['matched_class'].count()))
print('Number of annotations reweighted at the order level: ' +
      str(dt_isdb_results['matched_order'].count()))
print('Number of annotations reweighted at the family level: ' +
      str(dt_isdb_results['matched_family'].count()))
print('Number of annotations reweighted at the genus level: ' +
      str(dt_isdb_results['matched_genus'].count()))
print('Number of annotations reweighted at the species level: ' +
      str(dt_isdb_results['matched_species'].count()))


# %%

# we set the spectral score column as float
dt_isdb_results["score_input"] = pd.to_numeric(
    dt_isdb_results["score_input"], downcast="float")
# and we add it to the max txo score :
dt_isdb_results['score_input_taxo'] = dt_isdb_results['score_taxo'] + \
    dt_isdb_results['score_input']

# %%

dt_isdb_results['rank_spec_taxo'] = dt_isdb_results.groupby(
    'feature_id')['score_input_taxo'].rank(method='dense', ascending=False)

dt_isdb_results = dt_isdb_results.groupby(["feature_id"]).apply(
    lambda x: x.sort_values(["rank_spec_taxo"], ascending=True)).reset_index(drop=True)

# %%
# Get cluster Chemical class

for col in ['Superclass_cf_DNP', 'Class_cf_DNP', 'Subclass_cf_DNP', 'Parent_Level_1_cf_DNP']:

    df = dt_isdb_results.copy()
    df = df.drop_duplicates(subset=['feature_id', col])
    if use_post_taxo == True:
        df = df[df["component_id"] != -1]
        df = df[df.rank_spec_taxo <= top_N_chemical_consistency]
        df = df.groupby(
            ["component_id", col]
        ).agg({'feature_id': 'count',
               'rank_spec_taxo': 'mean'}
              ).reset_index(
        ).rename(columns={'feature_id': (col + '_count'),
                          'rank_spec_taxo': ('rank_' + col + '_mean')}
                 ).merge(cluster_count, on='component_id', how='left')
    else:
        df = df[df.component_id != -1]
        df = df[df.rank_spec <= top_N_chemical_consistency]
        df = df.groupby(
            ["component_id", col]
        ).agg({'feature_id': 'count',
               'rank_spec': 'mean'}
              ).reset_index(
        ).rename(columns={'feature_id': (col + '_count'),
                          'rank_spec': ('rank_' + col + '_mean')}
                 ).merge(cluster_count, on='component_id', how='left')

    df[('freq_' + col)] = df[(col + '_count')] / df['ci_count']
    df[(col + '_score')] = df[('freq_' + col)] / \
        (df[('rank_' + col + '_mean')]**(0.5))
    df = df.sort_values(
        (col + '_score'), ascending=False
    ).drop_duplicates(['component_id']
                      ).rename(columns={col: (col + '_consensus')})
    dt_isdb_results = dt_isdb_results.merge(
        df[[(col + '_consensus'), ('freq_' + col), 'component_id']], on='component_id', how='left')

# %% Chemical consistency reweighting

dt_isdb_results['Superclass_cf_score'] = dt_isdb_results.apply(
    lambda x: 1 if x.Superclass_cf_DNP == x.Superclass_cf_DNP_consensus else 0, axis=1)
dt_isdb_results['Class_cf_score'] = dt_isdb_results.apply(
    lambda x: 2 if x.Class_cf_DNP == x.Class_cf_DNP_consensus else 0, axis=1)
dt_isdb_results['Subclass_cf_score'] = dt_isdb_results.apply(
    lambda x: 3 if x.Subclass_cf_DNP == x.Subclass_cf_DNP_consensus else 0, axis=1)
dt_isdb_results['Parent_Level_1_cf_score'] = dt_isdb_results.apply(
    lambda x: 4 if x.Parent_Level_1_cf_DNP == x.Parent_Level_1_cf_DNP_consensus else 0, axis=1)

dt_isdb_results['score_max_consistency'] = dt_isdb_results[[
    "Superclass_cf_score",
    "Class_cf_score",
    "Subclass_cf_score",
    "Parent_Level_1_cf_score",
]].max(axis=1)

dt_isdb_results['Final_score'] = dt_isdb_results['score_input'] + \
    dt_isdb_results['score_taxo'] + dt_isdb_results['score_max_consistency']

dt_isdb_results['rank_final'] = dt_isdb_results.groupby(
    'feature_id')['Final_score'].rank(method='dense', ascending=False)


# %%

print('Number of annotations reweighted at the superclass level: ' +
      str(len(dt_isdb_results[(dt_isdb_results['Superclass_cf_score'] == 1)])))
print('Number of annotations reweighted at the class level: ' +
      str(len(dt_isdb_results[(dt_isdb_results['Class_cf_score'] == 2)])))
print('Number of annotations reweighted at the subclass level: ' +
      str(len(dt_isdb_results[(dt_isdb_results['Subclass_cf_score'] == 3)])))
print('Number of annotations reweighted at the parent level: ' +
      str(len(dt_isdb_results[(dt_isdb_results['Parent_Level_1_cf_score'] == 4)])))

# %%


dt_isdb_results_chem_rew = dt_isdb_results.loc[(
    dt_isdb_results.rank_final <= int(top_to_output))]
dt_isdb_results_chem_rew[["feature_id", "rank_final", "component_id"]] = dt_isdb_results_chem_rew[[
    "feature_id", "rank_final", "component_id"]].apply(pd.to_numeric, downcast='signed', axis=1)
dt_isdb_results_chem_rew = dt_isdb_results_chem_rew.sort_values(
    ["feature_id", "rank_final"], ascending=(False, True))
dt_isdb_results_chem_rew = dt_isdb_results_chem_rew.astype(str)

# %%
# Here we would like to filter results when short IK are repeated for the same feature_id at the same final rank
# see issue (https://gitlab.com/tima5/taxoscorer/-/issues/23)
# used 
# dt_isdb_results_chem_rew = dt_isdb_results_chem_rew.drop_duplicates(subset=['feature_id', 'short_inchikey', 'rank_final'], keep='first')

dt_isdb_results_chem_rew = dt_isdb_results_chem_rew.drop_duplicates(subset=['feature_id', 'short_inchikey'], keep='first')

dt_isdb_results_chem_rew = dt_isdb_results_chem_rew.astype({'feature_id' : 'float'})
dt_isdb_results_chem_rew = dt_isdb_results_chem_rew.astype({'feature_id' : 'int64'})


# %%

annot_attr = ['rank_spec', 'score_input', 'inchikey', 'libname', 'InChI_DNP',
              'Molecule_Name_DNP', 'Molecule_Formula_DNP', 'Accurate_Mass_DNP', 'matched_kingdom', 'matched_phylum', 'matched_class', 'matched_order',
              'matched_family', 'matched_genus', 'matched_species', 'score_taxo', 'score_max_consistency', 'Final_score', 'rank_final']

comp_attr = ['component_id', 'Superclass_cf_DNP_consensus', 'freq_Superclass_cf_DNP', 'Class_cf_DNP_consensus',
             'freq_Class_cf_DNP', 'Subclass_cf_DNP_consensus', 'freq_Subclass_cf_DNP', 'Parent_Level_1_cf_DNP_consensus',
             'freq_Parent_Level_1_cf_DNP']

col_to_keep = ['feature_id'] + comp_attr + annot_attr

df4cyto = dt_isdb_results_chem_rew[col_to_keep]

# %%

gb_spec = {c: '|'.join for c in annot_attr}
for c in comp_attr:
    gb_spec[c] = 'first'

# %%

df4cyto = df4cyto.groupby('feature_id').agg(gb_spec)

# %%
df4cyto.to_csv(output_weighed_ISDB_path, sep='\t')


# %%

print('Finished in %s seconds.' % (time.time() - start_time))
print('You can check your results here %s' % output_weighed_ISDB_path)


# %%

# dt_isdb_results = dt_isdb_results.astype(str)
# df4cyto = dt_isdb_results.groupby(['feature_id']).agg('|'.join)
# df4cyto.reset_index(drop=True, inplace=True)


# %%

# Testing cytoscape ready formatting



# df4cyto['rank_spec'] = df4cyto['rank_spec'].apply(lambda x: [x])

