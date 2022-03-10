import numpy as np
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import itertools

def plotter_single(dt_isdb_results_int, dt_samples_metadata,organism_header,treemap_chemo_counted_results_path, treemap_chemo_intensity_results_path):

    """ This function will get the CHEMBL ids from the structure_inchikey field

    Args:
        dt_isdb_results (dataframe) : a annotation table
        top_to_output (integer): Top N of candidate to keep
    Returns:
        dt_isdb_results_chem_rew (dataframe): a dataframe with the top N annotation ordered by final rank
    """

    dt_isdb_results_int = dt_isdb_results_int.replace({np.nan:'None'})


    unique_group_labels = dt_samples_metadata[organism_header].unique()

    pattern=[{"type": "domain"}]

    rep_pattern = list(itertools.chain.from_iterable(itertools.repeat(x, len(unique_group_labels)) for x in pattern))




    fig = make_subplots(1, len(unique_group_labels),
    subplot_titles = (unique_group_labels),
    specs=[rep_pattern])

    i=1
    for n in unique_group_labels:
        print(n)

        dt = dt_isdb_results_int[dt_isdb_results_int[n] > 0]
        fig.add_trace(px.treemap(dt, path=[px.Constant("all"), 'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class'], 
        values='counter').data[0], 
        row=1,col=i)
        i+=1

    # fig.update_traces(root_color="lightgrey")
    fig.update_layout(margin = dict(t=50, l=25, r=25, b=25),
    title_text="Metabolite annotation overview (size proportional to individual count)")
    fig.update_annotations(font_size=12)
    fig.show()
    fig.write_html(treemap_chemo_counted_results_path,
                full_html=False,
                include_plotlyjs='cdn')

    fig = make_subplots(1, len(unique_group_labels),
    subplot_titles = (unique_group_labels),
    specs=[rep_pattern])

    i=1
    for n in unique_group_labels:
        print(n)

        dt = dt_isdb_results_int[dt_isdb_results_int[n] > 0]
        fig.add_trace(px.treemap(dt, path=[px.Constant("all"), 'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class'], 
        values=n).data[0], 
        row=1,col=i)
        i+=1

    # fig.update_traces(root_color="lightgrey")
    fig.update_layout(margin = dict(t=50, l=25, r=25, b=25),
    title_text="Metabolite annotation overview (size proportional to mean intensity)")
    fig.update_annotations(font_size=12)
    fig.show()
    fig.write_html(treemap_chemo_intensity_results_path,
                full_html=False,
                include_plotlyjs='cdn')


def plotter_multi(dt_isdb_results_int, dt_samples_metadata, organism_header, sampletype_header, treemap_chemo_multi_counted_results_path, treemap_chemo_multi_intensity_results_path):

    """ This function will get the CHEMBL ids from the structure_inchikey field

    Args:
        dt_isdb_results (dataframe) : a annotation table
        top_to_output (integer): Top N of candidate to keep
    Returns:
        dt_isdb_results_chem_rew (dataframe): a dataframe with the top N annotation ordered by final rank
    """
    dt_isdb_results_int = dt_isdb_results_int.replace({np.nan:'None'})


    dt_samples_metadata['combined'] = dt_samples_metadata[organism_header] + '_' + dt_samples_metadata[sampletype_header]
    unique_group_labels = dt_samples_metadata['combined'].unique()
    type(unique_group_labels)

    pattern=[{"type": "domain"}]


    rep_pattern = list(itertools.chain.from_iterable(itertools.repeat(x, len(unique_group_labels)) for x in pattern))

    fig = make_subplots(1, len(unique_group_labels),
    subplot_titles = (unique_group_labels),
    specs=[rep_pattern])

    i=1
    for n in unique_group_labels:
        #print(n)

        dt = dt_isdb_results_int[dt_isdb_results_int[n] > 0]
        fig.add_trace(px.treemap(dt, path=[px.Constant("all"), 'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class'], 
        values='counter').data[0], 
        row=1,col=i)
        i+=1

    # fig.update_traces(root_color="lightgrey")
    fig.update_layout(margin = dict(t=50, l=25, r=25, b=25),
    title_text="Metabolite annotation overview (size proportional to individual count)")
    fig.update_annotations(font_size=12)
    fig.show()
    fig.write_html(treemap_chemo_multi_counted_results_path,
                full_html=False,
                include_plotlyjs='cdn')


    fig = make_subplots(1, len(unique_group_labels),
    subplot_titles = (unique_group_labels),
    specs=[rep_pattern])

    i=1
    for n in unique_group_labels:
        #print(n)

        dt = dt_isdb_results_int[dt_isdb_results_int[n] > 0]
        fig.add_trace(px.treemap(dt, path=[px.Constant("all"), 'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class'], 
        values=n).data[0], 
        row=1,col=i)
        i+=1

    # fig.update_traces(root_color="lightgrey")
    fig.update_layout(margin = dict(t=50, l=25, r=25, b=25),
    title_text="Metabolite annotation overview (size proportional to mean intensity)")
    fig.update_annotations(font_size=12)
    fig.show()
    fig.write_html(treemap_chemo_multi_intensity_results_path,
                full_html=False,
                include_plotlyjs='cdn')