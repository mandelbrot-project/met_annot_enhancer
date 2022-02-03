import pandas as pd
from tqdm import tqdm
from tqdm import tqdm_notebook


# We deactivate the iloc warning see https://stackoverflow.com/a/20627316
pd.options.mode.chained_assignment = None  # default='warn'


def ms1_matcher(input_df, adducts_file_path, ppm_tol, df_metadata):

    """Proceeds to MS1 match on an input dataframe using the precalculated adducts_dataframes and the appends metadata information
    Args:
        input_df (string) : the yaml parameters
        adducts_file_path (): 
        df_metadata ():
    Returns:
        df_MS1_merge (dataframe) :  a df with ms1 match according to a given tolerance and their annotations with an external metadat file
    """

    # Import MS1 list

    print('''
    Importing MS1 adducts file
    ''')

    adducts_df = pd.read_csv(
        adducts_file_path, compression='gzip', sep='\t')

    adducts_df['min'] = adducts_df['adduct_mass'] - \
        int(ppm_tol) * (adducts_df['adduct_mass'] / 1000000)
    adducts_df['max'] = adducts_df['adduct_mass'] + \
        int(ppm_tol) * (adducts_df['adduct_mass'] / 1000000)

    super_df = []

    for i, row in tqdm(input_df.iterrows(), total=input_df.shape[0]):
        par_mass = input_df.loc[i, 'mz']
        df_0 = input_df.loc[[i], [
            'feature_id', 'mz', 'component_id']]
        df_1 = adducts_df[(adducts_df['min'] <= par_mass) &
                          (adducts_df['max'] >= par_mass)]
        df_1['key'] = i
        df_1.drop(['min', 'max'], axis=1, inplace=True)
        df_tot = pd.merge(df_0, df_1, left_index=True,
                          right_on='key', how='left')
        super_df.append(df_tot)

    df_MS1 = pd.concat(super_df, axis=0)
    del super_df

    df_MS1 = df_MS1.drop(['key'], axis=1).drop_duplicates(
        subset=['feature_id', 'adduct'])

    df_MS1['libname'] = 'MS1_match'

    print('''
    MS1 annotation done !
    ''')

    df_meta_short = df_metadata[['short_inchikey', 'structure_exact_mass']]
    df_meta_short = df_meta_short.dropna(subset=['structure_exact_mass'])
    df_meta_short = df_meta_short.drop_duplicates(
        subset=['short_inchikey', 'structure_exact_mass'])

    df_meta_short = df_meta_short.round({'structure_exact_mass': 5})
    df_MS1 = df_MS1.round({'exact_mass': 5})

    df_MS1_merge = pd.merge(df_MS1, df_meta_short, left_on='exact_mass',
                            right_on='structure_exact_mass', how='left')
    df_MS1_merge = df_MS1_merge.dropna(subset=['short_inchikey'])

    df_MS1_merge['match_mzerror_MS1'] = df_MS1_merge['mz'] - \
        df_MS1_merge['adduct_mass']
    df_MS1_merge = df_MS1_merge.round({'match_mzerror_MS1': 5}).astype({
        'match_mzerror_MS1': 'str'})

    df_MS1_merge = df_MS1_merge.drop(
        ['structure_exact_mass', 'adduct_mass', 'exact_mass'], axis=1)
    df_MS1_merge['msms_score'] = 0

    return df_MS1_merge
