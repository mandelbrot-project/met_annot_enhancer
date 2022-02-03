# A set of loading functions
import pandas as pd 
import os

def isdb_results_loader(isdb_results_path):

    dt_isdb_results = pd.read_csv(isdb_results_path,
                                  sep='\t',
                                  usecols=['msms_score', 'feature_id',
                                           'reference_id', 'short_inchikey'],
                                  error_bad_lines=False, low_memory=True)

    ## we add a fixed in silico libname (to be changed later on)

    dt_isdb_results['libname'] = 'ISDB'

    return dt_isdb_results


def clusterinfo_summary_loader(clusterinfo_summary_path):

    ## we get the networks info (cluster id, component index and parent mass form the downloaded folder)

    clusterinfo_summary = pd.read_csv(clusterinfo_summary_path + str(os.listdir(clusterinfo_summary_path)[0]),
                                      sep='\t',
                                      usecols=['cluster index',
                                               'componentindex', 'parent mass'],
                                      error_bad_lines=False, low_memory=True)

    clusterinfo_summary.rename(columns={'cluster index': 'feature_id', 'componentindex': 'component_id',
                                        'parent mass': 'mz'}, inplace=True)

    # cluster_count = clusterinfo_summary.drop_duplicates(
    #     subset=['feature_id', 'component_id']).groupby("component_id").count()
    # cluster_count = cluster_count[['feature_id']].rename(
    #     columns={'feature_id': 'ci_count'}).reset_index()

    return clusterinfo_summary


def isdb_metadata_loader(isdb_metadata_path):

    df_metadata = pd.read_csv(isdb_metadata_path,
                            sep=',', error_bad_lines=False, low_memory=False)

    df_metadata['short_inchikey'] = df_metadata.structure_inchikey.str.split(
        "-", expand=True)[0]
    df_metadata.reset_index(inplace=True)

    # at this step we can only keep unique short_ik - organisms pairs

    df_metadata.drop_duplicates(subset=['short_inchikey', 'organism_wikidata'], keep='first', inplace=True, ignore_index=True)

    return df_metadata
