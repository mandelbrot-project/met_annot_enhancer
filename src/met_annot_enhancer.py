# %% Load packages

import pandas as pd
import numpy as np
import zipfile
import glob
import os
import shlex
import subprocess
from tqdm import tqdm
from tqdm import tqdm_notebook


# %% defininbg display options

pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 100)


# %% Choose line by line or classic (1 species) reweighting

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

# GNPS FBMN Job ID
job_id = "56d01c6ccfe143eca5252017202c8fef"
# toy set : f840934ade8744f9a6b4438804af8287
# LP: 75c64084d1e749748d5a29f671944231

# this is the path for the ISDB resukts file (topX)
#isdb_results_path = "J:/COMMON FASIE-FATHO/PF_project/Toy_Set/Toy_set_top_50_pos.out"
isdb_results_path = "/Users/pma/Dropbox/Research_UNIGE/git_repos/taxoscorer/data_in/hysope/hysope_pos_matchmsed_ISDB_DNP.tsv"
# isdb_results_path = "/Users/pma/tmp/bafu_ecometabo/FBMN_bafu_ecometabo_pos/FBMN_bafu_ecometabo_pos_msmatched_ISDB_DNP.out"


# this is the path for the DNP or OPENNPDB datatable file
metadata_path = "/Users/pma/Documents/190602_DNP_TAXcof_CF.tsv"

#Path to weighed annotation result of ISDB
output_weighed_ISDB_path = "/Users/pma/Dropbox/Research_UNIGE/git_repos/taxoscorer/data_in/hysope/hysope_pos_matchmsed_ISDB_DNP_repond.tsv"


# Path for the GNPS job export dir 

gnps_job_path = "/Users/pma/tmp/Fred_Legendre/"


# Set True if you want to use rank after taxonomical reweighting for consensus chemical class determination
use_post_taxo = True

# MS filename extension (a common pattern in all your filenames)
file_extension = '.mzXML'

# Set True if you want to use rank after taxonomical reweighting for consensus chemical class determination
top_N_chemical_consistency = 30

top_to_output = 3

ppm_tol = 5

#polarity = 'Neg'
polarity = 'Pos'

# %% Import MS1 list

if polarity == 'Pos':
    adducts_df = pd.read_csv(
        '../data_loc/db_pos.tsv.gz', compression='gzip', sep='\t')
else:
    adducts_df = pd.read_csv(
        '../data_loc/db_neg.tsv.gz', compression='gzip', sep='\t')

adducts_df['min'] = adducts_df['adduct_mass'] - \
    ppm_tol * (adducts_df['adduct_mass'] / 1000000)
adducts_df['max'] = adducts_df['adduct_mass'] + \
    ppm_tol * (adducts_df['adduct_mass'] / 1000000)


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

clusterinfo_summary = pd.read_csv('../../data_in/hysope/clusterinfo_summary/' + str(os.listdir('../../data_in/hysope/clusterinfo_summary/')[0]),
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



# %% Extracting biosource / feature for line by line
if Run_line_x_line == True:

    feature_intensity = pd.read_csv('../../data_in/hysope/quantification_table_reformatted/' + str(
        os.listdir('../../data_in/hysope/quantification_table_reformatted/')[0]), sep=',')
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

    Samples_metadata = pd.read_csv('../../data_in/hysope/metadata_table/' + str(os.listdir('../../data_in/hysope/metadata_table/')[0]), sep='\t',
                                   # usecols=['filename','ATTRIBUTE_phylum_cof', 'ATTRIBUTE_kingdom_cof',  'ATTRIBUTE_class_cof', 'ATTRIBUTE_order_cof', 'ATTRIBUTE_family_cof', 'ATTRIBUTE_genus_cof', 'ATTRIBUTE_species_cof'])
                                   usecols=['filename', 'ATTRIBUTE_Phylum', 'ATTRIBUTE_Kingdom',  'ATTRIBUTE_Class', 'ATTRIBUTE_Order', 'ATTRIBUTE_Family', 'ATTRIBUTE_Genus', 'ATTRIBUTE_Species'])

    res = feature_intensity[feature_intensity != 0].stack()
    df_res = res.to_frame().reset_index()
    df_merged = pd.merge(df_res, Samples_metadata, left_on='MS_filename',
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
       dt_isdb_results['ATTRIBUTE_species_cof'] = species_bio
       dt_isdb_results['ATTRIBUTE_genus_cof'] = genus_bio
       dt_isdb_results['ATTRIBUTE_family_cof'] = family_bio
       dt_isdb_results['ATTRIBUTE_order_cof'] = order_bio
       dt_isdb_results['ATTRIBUTE_class_cof'] = class_bio
       dt_isdb_results['ATTRIBUTE_phylum_cof'] = phylum_bio
       dt_isdb_results['ATTRIBUTE_kingdom_cof'] = kingdom_bio
       
#%% Taxonomical Reweighting

cols_ref = ['Kingdom_cof_DNP', 'Phylum_cof_DNP',  'Class_cof_DNP', 'Order_cof_DNP', 'Family_cof_DNP', 'Genus_cof_DNP', 'Species_cof_DNP']
cols_att = ['ATTRIBUTE_Kingdom', 'ATTRIBUTE_Phylum',  'ATTRIBUTE_Class', 'ATTRIBUTE_Order', 'ATTRIBUTE_Family', 'ATTRIBUTE_Genus', 'ATTRIBUTE_Species']
#cols_att = ['ATTRIBUTE_kingdom_cof', 'ATTRIBUTE_phylum_cof', 'ATTRIBUTE_class_cof', 'ATTRIBUTE_order_cof', 'ATTRIBUTE_family_cof', 'ATTRIBUTE_genus_cof', 'ATTRIBUTE_species_cof']
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
    dt_isdb_results.rank_final <= top_to_output)]
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

# dt_isdb_results = dt_isdb_results.astype(str)
# df4cyto = dt_isdb_results.groupby(['feature_id']).agg('|'.join)
# df4cyto.reset_index(drop=True, inplace=True)


# %%

# Testing cytoscape ready formatting



# df4cyto['rank_spec'] = df4cyto['rank_spec'].apply(lambda x: [x])

