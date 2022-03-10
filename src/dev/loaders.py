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


def isdb_metadata_loader(isdb_metadata_path, organism_header):

    df_metadata = pd.read_csv(isdb_metadata_path,
                              sep=',', error_bad_lines=False, low_memory=False)

    df_metadata['short_inchikey'] = df_metadata.structure_inchikey.str.split(
        "-", expand=True)[0]
    df_metadata.reset_index(inplace=True)

    # at this step we can only keep unique short_ik - organisms pairs

    df_metadata.drop_duplicates(subset=[
                                'short_inchikey', organism_header], keep='first', inplace=True, ignore_index=True)

    return df_metadata


def samples_metadata_loader(samples_metadata_table_path, organism_header):
    # the metadata table is loaded using the organism column specified before
    samples_metadata = pd.read_csv(samples_metadata_table_path + str(os.listdir(samples_metadata_table_path)[0]), sep='\t',
                                   usecols=['filename', organism_header])
    return samples_metadata

def samples_metadata_full_loader(samples_metadata_table_path):
    # the metadata table is loaded using the organism column specified before
    samples_metadata = pd.read_csv(samples_metadata_table_path + str(os.listdir(samples_metadata_table_path)[0]), sep='\t')
    
    return samples_metadata


def feature_intensity_table_loader(feature_intensity_table_path):
    feature_intensity = pd.read_csv(feature_intensity_table_path + str(
        os.listdir(feature_intensity_table_path)[0]), sep=',')

    return feature_intensity
