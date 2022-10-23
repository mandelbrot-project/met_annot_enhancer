import pandas as pd
from loaders import feature_intensity_table_loader

def feature_intensity_table_formatter(feature_intensity_table_path, file_extension, msfile_suffix):

    """Formats a feature intensity table to an appropriate format for biosource_contribution_fetcher()
    Args:
        feature_intensity_table_path (str): a path to a feature intensity table 
        file_extension (str): a string to match the filename extension(typically .mzXML or similar) defines in the .yml file
        msfile_suffix (str): a string to match an eventual filename suffix (example Peak height)
    Returns:
        feature_intensity_table (dataframe): a formatted feature intensity table 
    """
    # The feature_intensity_table is loaded
    feature_intensity_table = feature_intensity_table_loader(feature_intensity_table_path)

    # formatting the feature table 
    feature_intensity_table.rename(columns={'row ID': 'row_ID'}, inplace=True)
    feature_intensity_table.set_index('row_ID', inplace=True)
    feature_intensity_table = feature_intensity_table.filter(regex=file_extension)
    feature_intensity_table.columns = feature_intensity_table.columns.str.replace(msfile_suffix, '') 
    # this above is not safe, we should find an alternative. Maybe raising an issue if the suffix is not found 
    feature_intensity_table.rename_axis("MS_filename", axis="columns", inplace = True)

    return feature_intensity_table


def cytoscape_attributes_formatter(df_input):

    """Formats a feature intensity table to an appropriate format for biosource_contribution_fetcher()
    Args:
        feature_intensity_table_path (str): a path to a feature intensity table 
        file_extension (str): a string to match the filename extension(typically .mzXML or similar) defines in the .yml file
        msfile_suffix (str): a string to match an eventual filename suffix (example Peak height)
    Returns:
        feature_intensity_table (dataframe): a formatted feature intensity table 
    """

    all_columns = list(df_input) # Creates list of all column headers
    
    df_input[all_columns] = df_input[all_columns].astype(str)

    gb_spec = {c: '|'.join for c in annot_attr}

    for c in comp_attr:
        gb_spec[c] = 'first'

    df_output = df_input.groupby('feature_id').agg(gb_spec)

    return df_output


def table_for_plots_formatter(df_flat, feature_intensity_table_formatted, dt_samples_metadata, organism_header, var_one_header, multi_plot):

    """Formats a feature intensity table to an appropriate format for biosource_contribution_fetcher()
    Args:
        feature_intensity_table_path (str): a path to a feature intensity table 
        file_extension (str): a string to match the filename extension(typically .mzXML or similar) defines in the .yml file
        msfile_suffix (str): a string to match an eventual filename suffix (example Peak height)
    Returns:
        feature_intensity_table (dataframe): a formatted feature intensity table 
    """


    feature_intensity_table_t = feature_intensity_table_formatted.transpose()

    feature_intensity_meta = pd.merge(left=dt_samples_metadata, right=feature_intensity_table_t, left_on='filename', right_on='MS_filename',how='inner')


    feature_intensity_meta_gp_species = feature_intensity_meta.groupby(organism_header).mean()
    feature_intensity_meta_gp_species = feature_intensity_meta_gp_species.transpose()
    feature_intensity_meta_gp_species.index.name = 'row_ID'


    feature_intensity_table_formatted.reset_index(inplace=True)
    feature_intensity_meta_gp_species.reset_index(inplace=True)


    ft_merged = pd.merge(feature_intensity_table_formatted, feature_intensity_meta_gp_species, on='row_ID', how='left')


    if multi_plot == True:
        # a security when numeric values are passed
        feature_intensity_meta_gp_multi = feature_intensity_meta.groupby([organism_header,var_one_header]).mean()
        feature_intensity_meta_gp_multi = feature_intensity_meta_gp_multi.transpose()
        feature_intensity_meta_gp_multi.columns = feature_intensity_meta_gp_multi.columns.map('_'.join)
        feature_intensity_meta_gp_multi.index.name = 'row_ID'
        feature_intensity_meta_gp_multi.reset_index(inplace=True)

        ft_merged = pd.merge(ft_merged, feature_intensity_meta_gp_multi, on='row_ID', how='left')


    # df_flat['feature_id'] = df_flat['feature_id'].astype('int')

    dt_isdb_results_int = pd.merge(
        df_flat, ft_merged, left_on='feature_id', right_on='row_ID', how='left')

    dt_isdb_results_int['counter'] = 1


    return dt_isdb_results_int




def samples_metadata_filterer(dt_samples_metadata, organism_header, var_one_header, drop_pattern):

    """Formats a feature intensity table to an appropriate format for biosource_contribution_fetcher()
    Args:
        feature_intensity_table_path (str): a path to a feature intensity table 
        file_extension (str): a string to match the filename extension(typically .mzXML or similar) defines in the .yml file
        msfile_suffix (str): a string to match an eventual filename suffix (example Peak height)
    Returns:
        feature_intensity_table (dataframe): a formatted feature intensity table 
    """

    if len(drop_pattern) != 0:
        dt_samples_metadata = dt_samples_metadata[~dt_samples_metadata[organism_header].str.contains(drop_pattern, na=False)]
        dt_samples_metadata = dt_samples_metadata[~dt_samples_metadata[var_one_header].str.contains(drop_pattern, na=False)]
    else:
        dt_samples_metadata = dt_samples_metadata

    return dt_samples_metadata


def samples_metadata_filterer_sampletype(dt_samples_metadata, organism_header, var_one_header, sampletype_header, sampletype_value_sample, drop_pattern, multi_plot):

    """Formats a feature intensity table to an appropriate format for biosource_contribution_fetcher()
    Args:
        feature_intensity_table_path (str): a path to a feature intensity table 
        file_extension (str): a string to match the filename extension(typically .mzXML or similar) defines in the .yml file
        msfile_suffix (str): a string to match an eventual filename suffix (example Peak height)
    Returns:
        feature_intensity_table (dataframe): a formatted feature intensity table 
    """

    if len(drop_pattern) != 0 and multi_plot == True:
        dt_samples_metadata = dt_samples_metadata[~dt_samples_metadata[organism_header].str.contains(drop_pattern, na=False)]
        dt_samples_metadata = dt_samples_metadata[~dt_samples_metadata[var_one_header].str.contains(drop_pattern, na=False)]
        dt_samples_metadata = dt_samples_metadata[dt_samples_metadata[sampletype_header] == sampletype_value_sample]
    if len(drop_pattern) != 0 and multi_plot == False:
        dt_samples_metadata = dt_samples_metadata[~dt_samples_metadata[organism_header].str.contains(drop_pattern, na=False)]
        dt_samples_metadata = dt_samples_metadata[dt_samples_metadata[sampletype_header] == sampletype_value_sample]
    else:
        dt_samples_metadata = dt_samples_metadata

    return dt_samples_metadata

