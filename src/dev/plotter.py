import numpy as np
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import itertools

# Microshades color palettes used for chemical classes
# r$> hex_custom
#   micro_cvd_gray micro_cvd_purple micro_cvd_blue micro_cvd_orange micro_cvd_green micro_cvd_turquoise micro_orange micro_purple
# 1        #616161          #7D3560        #098BD9          #9D654C         #4E7705             #148F77      #ff7f00      #6a51a3
# 2        #8B8B8B          #A1527F        #56B4E9          #C17754         #6D9F06             #009E73      #fe9929      #807dba
# 3        #B7B7B7          #CC79A7        #7DCCFF          #F09163         #97CE2F             #43BA8F      #fdae6b      #9e9ac8
# 4        #D6D6D6          #E794C1        #BCE1FF          #FCB076         #BDEC6F             #48C9B0      #fec44f      #bcbddc
# 5        #F5F5F5          #EFB6D6        #E7F4FF          #FFD5AF         #DDFFA0             #A3E4D7      #feeda0      #dadaeb

    #   "Terpenoids" = "micro_cvd_purple",
    #   "Fatty acids" = "micro_cvd_blue",
    #   "Polyketides" = "micro_cvd_orange",
    #   "Alkaloids" = "micro_cvd_green",
    #   "Shikimates and Phenylpropanoids" = "micro_cvd_turquoise",
    #   "Amino acids and Peptides" = "micro_orange",
    #   "Carbohydrates" = "micro_purple",
    #   "Other" = "micro_cvd_gray"


def plotter_single(dt_isdb_results_int, dt_samples_metadata,organism_header,treemap_chemo_counted_results_path, treemap_chemo_intensity_results_path):
    """
    Generate a treemap visualization for metabolite annotation data for a single variable.

    Args:
        dt_isdb_results_int (DataFrame): Annotation table with metabolite data.
        dt_samples_metadata (DataFrame): Metadata table with sample information.
        organism_header (str): Column name in dt_samples_metadata representing the organism labels.
        treemap_chemo_counted_results_path (str): File path to save the treemap visualization for counted results.
        treemap_chemo_intensity_results_path (str): File path to save the treemap visualization for intensity results.

    Returns:
        None
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
                     color='structure_taxonomy_npclassifier_01pathway',
                     color_discrete_map={
                         '(?)':'#F5F5F5',
                         'Terpenoids':'#7D3560',
                         'Alkaloids': '#4E7705',
                         'Amino acids and Peptides': '#ff7f00',
                         'Polyketides': '#9D654C',
                         'Shikimates and Phenylpropanoids': '#148F77',
                         'Fatty acids': '#098BD9',
                         'Carbohydrates': '#6a51a3',
                         'None': '#616161',},
        values='counter').data[0], 
        row=1,col=i)
        # fig.update_traces(marker_cornerradius=5, selector=dict(type='treemap'))
        i+=1

    # fig.update_traces(root_color="lightgrey")
    fig.update_layout(margin = dict(t=50, l=25, r=25, b=25),
    title_text="Metabolite annotation overview (size proportional to individual count)")
    fig.update_annotations(font_size=12)
    # fig.update_traces(marker=dict(cornerradius=5))
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
                     color='structure_taxonomy_npclassifier_01pathway',
                     color_discrete_map={
                         '(?)':'#F5F5F5',
                         'Terpenoids':'#7D3560',
                         'Alkaloids': '#4E7705',
                         'Amino acids and Peptides': '#ff7f00',
                         'Polyketides': '#9D654C',
                         'Shikimates and Phenylpropanoids': '#148F77',
                         'Fatty acids': '#098BD9',
                         'Carbohydrates': '#6a51a3',
                         'None': '#616161',},
        values=n).data[0], 
        row=1,col=i)
        i+=1

    # fig.update_traces(root_color="lightgrey")
    fig.update_layout(margin = dict(t=50, l=25, r=25, b=25),
    title_text="Metabolite annotation overview (size proportional to mean intensity)")
    fig.update_annotations(font_size=12)
    # fig.update_traces(marker=dict(cornerradius=5))
    fig.show()
    fig.write_html(treemap_chemo_intensity_results_path,
                full_html=False,
                include_plotlyjs='cdn')


def plotter_multi(dt_isdb_results_int, dt_samples_metadata, organism_header, var_one_header, treemap_chemo_multi_counted_results_path, treemap_chemo_multi_intensity_results_path):
    """
    Generate treemap visualizations for metabolite annotation data for multiple variables.

    Args:
        dt_isdb_results_int (DataFrame): Annotation table with metabolite data.
        dt_samples_metadata (DataFrame): Metadata table with sample information.
        organism_header (str): Column name in dt_samples_metadata representing the organism labels.
        var_one_header (str): Column name in dt_samples_metadata representing the first variable labels.
        treemap_chemo_multi_counted_results_path (str): File path to save the treemap visualization for counted results.
        treemap_chemo_multi_intensity_results_path (str): File path to save the treemap visualization for intensity results.

    Returns:
        None
    """
    dt_isdb_results_int = dt_isdb_results_int.replace({np.nan:'None'})


    dt_samples_metadata['combined'] = dt_samples_metadata[organism_header] + '_' + dt_samples_metadata[var_one_header]
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
                     color='structure_taxonomy_npclassifier_01pathway',
                     color_discrete_map={
                         '(?)':'#F5F5F5',
                         'Terpenoids':'#7D3560',
                         'Alkaloids': '#4E7705',
                         'Amino acids and Peptides': '#ff7f00',
                         'Polyketides': '#9D654C',
                         'Shikimates and Phenylpropanoids': '#148F77',
                         'Fatty acids': '#098BD9',
                         'Carbohydrates': '#6a51a3',
                         'None': '#616161',},
        values='counter').data[0], 
        row=1,col=i)
        i+=1

    # fig.update_traces(root_color="lightgrey")
    fig.update_layout(margin = dict(t=50, l=25, r=25, b=25),
    title_text="Metabolite annotation overview (size proportional to individual count)")
    fig.update_annotations(font_size=12)
    # fig.update_traces(marker=dict(cornerradius=5))
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
                     color='structure_taxonomy_npclassifier_01pathway',
                     color_discrete_map={
                         '(?)':'#F5F5F5',
                         'Terpenoids':'#7D3560',
                         'Alkaloids': '#4E7705',
                         'Amino acids and Peptides': '#ff7f00',
                         'Polyketides': '#9D654C',
                         'Shikimates and Phenylpropanoids': '#148F77',
                         'Fatty acids': '#098BD9',
                         'Carbohydrates': '#6a51a3',
                         'None': '#616161',},
        values=n).data[0], 
        row=1,col=i)
        fig.update_traces(root_color="whitesmoke")
        i+=1

    # fig.update_traces(root_color="lightgrey")
    fig.update_layout(margin = dict(t=50, l=25, r=25, b=25),
    title_text="Metabolite annotation overview (size proportional to mean intensity)")
    fig.update_annotations(font_size=12)
    # fig.update_traces(marker=dict(cornerradius=5))
    fig.show()
    fig.write_html(treemap_chemo_multi_intensity_results_path,
                full_html=False,
                include_plotlyjs='cdn')
