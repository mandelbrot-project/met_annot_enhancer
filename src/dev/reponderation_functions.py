import pandas as pd
import numpy as np

from helpers import cluster_counter

def biosource_contribution_fetcher(feature_intensity_table, samples_metadata, top_n):

    """Generates pathes used by the script according to parameters of the yaml file
    Args:
        params_list (list) : the yaml parameters
    Returns:
        pathes (str)
    """
    
    print('''
    Fetching the biosource contribution per feature ...
    ''')


    # fetching the TopN contributions
    if top_n == 0:
        feature_intensity_table = feature_intensity_table.where(feature_intensity_table.apply(
            lambda x: x.isin(x.nlargest(len(feature_intensity_table.columns))), axis=1), 0)  # top N here
    else:
        feature_intensity_table = feature_intensity_table.where(feature_intensity_table.apply(
            lambda x: x.isin(x.nlargest(top_n)), axis=1), 0)  # top N here


    res = feature_intensity_table[feature_intensity_table != 0].stack()
    df_res = res.to_frame().reset_index()


    # merging back with the samples metadata and aggregating as lists

    df_merged = pd.merge(df_res, samples_metadata, left_on='MS_filename',
                            right_on='filename', how='left').drop([0, 'MS_filename', 'filename'], axis=1)
    
    df_merged = df_merged[['row_ID','query_otol_domain', 'query_otol_kingdom', 'query_otol_phylum', 'query_otol_class',
            'query_otol_order', 'query_otol_family', 'query_otol_tribe', 'query_otol_genus', 'query_otol_species']]
    df_merged = df_merged.groupby('row_ID').agg(lambda x: list(x))
    df_merged.reset_index(inplace=True)

    return df_merged



def taxonomical_reponderator(dt_isdb_results, min_score_taxo_ms1):

    """Generates pathes used by the script according to parameters of the yaml file
    Args:
        params_list (list) : the yaml parameters
    Returns:
        pathes (str)
    """



    print('''
    Proceeding to taxonomically informed reponderation ...
    ''')

    cols_ref = ['organism_taxonomy_01domain', 'organism_taxonomy_02kingdom',  'organism_taxonomy_03phylum', 'organism_taxonomy_04class',
                'organism_taxonomy_05order', 'organism_taxonomy_06family', 'organism_taxonomy_07tribe', 'organism_taxonomy_08genus', 'organism_taxonomy_09species']

    cols_att = ['query_otol_domain', 'query_otol_kingdom', 'query_otol_phylum', 'query_otol_class',
                'query_otol_order', 'query_otol_family', 'query_otol_tribe', 'query_otol_genus', 'query_otol_species']

    cols_match = ['matched_domain', 'matched_kingdom', 'matched_phylum', 'matched_class',
                'matched_order', 'matched_family', 'matched_tribe', 'matched_genus', 'matched_species']

    col_prev = None

    for col_ref, col_att, col_match in zip(cols_ref, cols_att, cols_match):

            dt_isdb_results[col_ref].fillna('Unknown', inplace=True)
            dt_isdb_results[col_ref] = dt_isdb_results[col_ref].apply(lambda x: [x])
            dt_isdb_results[col_match] = [list(set(a).intersection(set(b))) for a, b in zip(dt_isdb_results[col_ref], dt_isdb_results[col_att])] # Allows to compare 2 lists
            dt_isdb_results[col_match] = dt_isdb_results[col_match].apply(lambda y: np.nan if len(y)==0 else y)

            if col_prev != None:
                    dt_isdb_results[col_match].where(dt_isdb_results[col_prev].notnull(), np.nan)
            col_prev = col_match

    # Note for future self. If you get a TypeError: unhashable type: 'list' error. before messing around with the previous line make sure that the taxonomy has been appended at the dt_isdb_results = pd.merge(
    #  ' dt_isdb_results, df_merged, left_on='feature_id', right_on='row_ID', how='left')' step before. Usuall this comes from a bad definition of the regex (ex .mzXMl insted of .mzML) in the params file. Should find a safer way to deal with these extensions in the header.


    dt_isdb_results['score_taxo'] = dt_isdb_results[cols_match].count(axis=1)

    # Filter out MS1 annotations without a reweighting at a given taxo level prior to chemo repond


    dt_isdb_results = dt_isdb_results[
        (dt_isdb_results['score_taxo'] >= min_score_taxo_ms1) | (
        dt_isdb_results['libname'] == 'ISDB')]


    # we set the spectral score column as float
    dt_isdb_results["msms_score"] = pd.to_numeric(
        dt_isdb_results["msms_score"], downcast="float")
    # and we add it to the max txo score :
    dt_isdb_results['msms_score_taxo'] = dt_isdb_results['score_taxo'] + \
        dt_isdb_results['msms_score']


    dt_isdb_results['rank_spec_taxo'] = dt_isdb_results.groupby(
        'feature_id')['msms_score_taxo'].rank(method='dense', ascending=False)

    dt_isdb_results = dt_isdb_results.groupby(["feature_id"]).apply(
        lambda x: x.sort_values(["rank_spec_taxo"], ascending=True)).reset_index(drop=True)


    print('Total number of annotations after filtering MS1 annotations not reweighted at the minimal taxonomical level: ' +
        str(len(dt_isdb_results)))

    print('Number of annotations reweighted at the domain level: ' +
        str(dt_isdb_results['matched_domain'].count()))
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
    print('Number of annotations reweighted at the tribe level: ' +
        str(dt_isdb_results['matched_tribe'].count()))
    print('Number of annotations reweighted at the genus level: ' +
        str(dt_isdb_results['matched_genus'].count()))
    print('Number of annotations reweighted at the species level: ' +
        str(dt_isdb_results['matched_species'].count()))

    return dt_isdb_results



def chemical_reponderator(clusterinfo_summary_file, dt_isdb_results, top_N_chemical_consistency):

    """Generates pathes used by the script according to parameters of the yaml file
    Args:
        params_list (list) : the yaml parameters
    Returns:
        pathes (str)

    """

    cluster_count = cluster_counter(clusterinfo_summary_file)


    # Get cluster Chemical class

    for col in ['structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class']:

        df = dt_isdb_results.copy()
        df = df.drop_duplicates(subset=['feature_id', col])
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

        df[('freq_' + col)] = df[(col + '_count')] / df['ci_count']
        df[(col + '_score')] = df[('freq_' + col)] / \
            (df[('rank_' + col + '_mean')]**(0.5))

        df = df.sort_values(
            (col + '_score'), ascending=False
        ).drop_duplicates(['component_id']
                          ).rename(columns={col: (col + '_consensus')})

        dt_isdb_results = dt_isdb_results.merge(
            df[[(col + '_consensus'), ('freq_' + col), 'component_id']], on='component_id', how='left')

    # Chemical consistency reweighting

    print('''
    Proceeding to chemically informed reponderation ...
    ''')

    dt_isdb_results['structure_taxonomy_npclassifier_01pathway_score'] = dt_isdb_results.apply(
        lambda x: 1 if x.structure_taxonomy_npclassifier_01pathway == x.structure_taxonomy_npclassifier_01pathway_consensus else 0, axis=1)
    dt_isdb_results['structure_taxonomy_npclassifier_02superclass_score'] = dt_isdb_results.apply(
        lambda x: 2 if x.structure_taxonomy_npclassifier_02superclass == x.structure_taxonomy_npclassifier_02superclass_consensus else 0, axis=1)
    dt_isdb_results['structure_taxonomy_npclassifier_03class_score'] = dt_isdb_results.apply(
        lambda x: 3 if x.structure_taxonomy_npclassifier_03class == x.structure_taxonomy_npclassifier_03class_consensus else 0, axis=1)

    dt_isdb_results['score_max_consistency'] = dt_isdb_results[[
        "structure_taxonomy_npclassifier_01pathway_score",
        "structure_taxonomy_npclassifier_02superclass_score",
        "structure_taxonomy_npclassifier_03class_score"
    ]].max(axis=1)

    dt_isdb_results['final_score'] = dt_isdb_results['msms_score'] + \
        dt_isdb_results['score_taxo'] + \
        dt_isdb_results['score_max_consistency']

    dt_isdb_results['rank_final'] = dt_isdb_results.groupby(
        'feature_id')['final_score'].rank(method='dense', ascending=False)

    print('Number of annotations reweighted at the NPClassifier pathway level: ' +
          str(len(dt_isdb_results[(dt_isdb_results['structure_taxonomy_npclassifier_01pathway_score'] == 1)])))
    print('Number of annotations reweighted at the NPClassifier superclass level: ' +
          str(len(dt_isdb_results[(dt_isdb_results['structure_taxonomy_npclassifier_02superclass_score'] == 2)])))
    print('Number of annotations reweighted at the NPClassifier class level: ' +
          str(len(dt_isdb_results[(dt_isdb_results['structure_taxonomy_npclassifier_03class_score'] == 3)])))

    return dt_isdb_results



