
import pandas as pd
from opentree import OT
import json
from pandas import json_normalize
import numpy as np


def tnrs_fetcher(species, path_to_results_folders, project_name):

    species_tnrs_matched = OT.tnrs_match(
        species, context_name=None, do_approximate_matching=True, include_suppressed=False)

    with open(str(path_to_results_folders + project_name + '_' + 'species.json'), 'w') as out:
        sf = json.dumps(species_tnrs_matched.response_dict,
                        indent=2, sort_keys=True)
        out.write('{}\n'.format(sf))


def taxon_info_fetcher(ott_list, path_to_results_folders, project_name):

    taxon_info = []

    for i in ott_list:
        query = OT.taxon_info(i, include_lineage=True)
        taxon_info.append(query)

    tl = []

    for i in taxon_info:
        with open(str(path_to_results_folders + project_name + '_' + 'taxon_info.json'), 'w') as out:
            tl.append(i.response_dict)
            yo = json.dumps(tl)
            out.write('{}\n'.format(yo))


def taxa_lineage_appender(samples_metadata, organism_header, do_taxo_resolving, path_to_results_folders, project_name):

    # Now we want to get the taxonomic information for each of the samples
    # so we first want to extract the species information from the metadata file

    samples_metadata[organism_header].dropna(inplace=True)
    samples_metadata[organism_header] = samples_metadata[organism_header].str.lower()
    species = samples_metadata[organism_header].unique()
    len_species = len(species)

    print("%s unique species have been selected from the metadata table." % len_species)

    if do_taxo_resolving == True:

        tnrs_fetcher(species, path_to_results_folders, project_name)

    with open(str(path_to_results_folders + project_name + '_' + 'species.json')) as tmpfile:
        jsondic = json.loads(tmpfile.read())

    json_normalize(jsondic)

    df_species_tnrs_matched = json_normalize(jsondic,
                                             record_path=['results', 'matches']
                                             )
    df_species_tnrs_unmatched = json_normalize(jsondic,
                                               record_path=['unmatched_names']
                                               )

    # We then want to match with the accepted name instead of the synonym in case both are present.
    # We thus order by matched_name and then by is_synonym status prior to returning the first row.

    df_species_tnrs_matched.sort_values(
        ['search_string', 'is_synonym'], axis=0, inplace=True)
    df_species_tnrs_matched_unique = df_species_tnrs_matched.drop_duplicates(
        'search_string', keep='first')

    # both df are finally merged
    merged_df = pd.merge(samples_metadata, df_species_tnrs_matched_unique, how='left',
                         left_on=organism_header, right_on='search_string', indicator=True)

    # converting 'ott_ids' from float to int (check the astype('Int64') whic will work while the astype('int') won't see https://stackoverflow.com/a/54194908)
    merged_df['taxon.ott_id'] = merged_df['taxon.ott_id'].astype('Int64')

    # However, we then need to put them back to
    merged_df['taxon.ott_id']
    ott_list = list(merged_df['taxon.ott_id'].dropna().astype('int'))

    if do_taxo_resolving == True:

        taxon_info_fetcher(ott_list, path_to_results_folders, project_name)

    with open(str(path_to_results_folders + project_name + '_' + 'taxon_info.json')) as tmpfile:
        jsondic = json.loads(tmpfile.read())

    df = json_normalize(jsondic)

    df_tax_lineage = json_normalize(jsondic,
                                    record_path=['lineage'],
                                    meta=['ott_id', 'unique_name'],
                                    record_prefix='sub_',
                                    errors='ignore'
                                    )

    # This keeps the last occurence of each ott_id / sub_rank grouping https://stackoverflow.com/a/41886945
    df_tax_lineage_filtered = df_tax_lineage.groupby(
        ['ott_id', 'sub_rank'], as_index=False).last()

    #Here we pivot long to wide to get the taxonomy
    df_tax_lineage_filtered_flat = df_tax_lineage_filtered.pivot(
        index='ott_id', columns='sub_rank', values='sub_name')

    # Here we actually also want the lowertaxon (species usually) name
    df_tax_lineage_filtered_flat = pd.merge(df_tax_lineage_filtered_flat, df_tax_lineage_filtered[[
                                            'ott_id', 'unique_name']], how='left', on='ott_id', )

    #Despite the left join ott_id are duplicated
    df_tax_lineage_filtered_flat.drop_duplicates(
        subset=['ott_id', 'unique_name'], inplace=True)

    # here we want to have these columns whatevere happens
    col_list = ['ott_id', 'domain', 'kingdom', 'phylum',
                'class', 'order', 'family', 'tribe', 'genus', 'unique_name']

    df_tax_lineage_filtered_flat = df_tax_lineage_filtered_flat.reindex(
        columns=col_list, fill_value=np.NaN)

    # We now rename our columns of interest
    renaming_dict = {'domain': 'query_otol_domain',
                     'kingdom': 'query_otol_kingdom',
                     'phylum': 'query_otol_phylum',
                     'class': 'query_otol_class',
                     'order': 'query_otol_order',
                     'family': 'query_otol_family',
                     'tribe': 'query_otol_tribe',
                     'genus': 'query_otol_genus',
                     'unique_name': 'query_otol_species'}

    df_tax_lineage_filtered_flat.rename(columns=renaming_dict, inplace=True)

    # We select columns of interest
    cols_to_keep = ['ott_id',
                    'query_otol_domain',
                    'query_otol_kingdom',
                    'query_otol_phylum',
                    'query_otol_class',
                    'query_otol_order',
                    'query_otol_family',
                    'query_otol_tribe',
                    'query_otol_genus',
                    'query_otol_species']

    df_tax_lineage_filtered_flat = df_tax_lineage_filtered_flat[cols_to_keep]

    # We merge this back with the samplemetadata only if we have an ott.id in the merged df
    samples_metadata = pd.merge(merged_df[pd.notnull(merged_df['taxon.ott_id'])],
                                df_tax_lineage_filtered_flat, how='left', left_on='taxon.ott_id', right_on='ott_id')

    return samples_metadata
