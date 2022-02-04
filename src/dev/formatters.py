
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


