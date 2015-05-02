import logging

import numpy

from matplotlib import pyplot

PIE_CHART_COLORS = [
    'yellowgreen',
    'gold',
    'lightcoral',
    'lightskyblue',
    'red',
    'orange',
    'purple',
    'brown',
    'black',
    'beige',
]

def abbreviate_allele(allele):
    if allele == "":
        return "del"
    if len(allele) > 5:
        return "%s+%d" % (allele[0], len(allele))
    return allele

def text_args(user_args, **default_args):
    if not isinstance(user_args, dict):
        user_args = {"s": user_args}
    default_args.update(user_args)
    return default_args

def plot(
        allele_vectors_matrix,
        radii,
        row_labels,
        col_labels,
        cell_labels,
        rows_per_figure=10,
        col_label_height=8.0,
        row_label_width=3.5,
        cell_width=2.5,
        cell_height=3.5,
        pie_chart_colors=PIE_CHART_COLORS,
        legend_extra_properties={},
        subplots_adjust_properties={'hspace': 0.1}):
    """

    Parameters
    ---

    matrix_allele_dicts
        List of DataFrame instances. Each dataframe has rows equal to
        number of columns in the plot. The columns are the alleles.

    radii
        Num rows x num cols array

    """

    num_rows = len(radii)
    num_cols = len(radii[0])

    assert len(allele_vectors_matrix) == num_rows
    assert len(row_labels) == num_rows
    assert len(col_labels) == num_cols

    legend_properties = {
        'bbox_to_anchor': (-.5, .8),
        'fontsize': 'large',
        'loc': 'upper left',
        'frameon': False,
    }
    legend_properties.update(legend_extra_properties)

    figure_row_plans = []
    current_figure_row_end = 0  # exclusive
    while current_figure_row_end < num_rows:
        current_figure_rows = min(
            rows_per_figure, (num_rows - current_figure_row_end))
        current_figure_row_start = current_figure_row_end
        current_figure_row_end += current_figure_rows
        figure_row_plans.append((
            current_figure_rows,
            current_figure_row_start,
            current_figure_row_end))

    for (figure_num, plan) in enumerate(figure_row_plans):
        logging.info("Creating figure %d / %d" % (
            figure_num + 1, len(figure_row_plans)))
        (current_figure_rows,
            current_figure_row_start,
            current_figure_row_end) = plan
   
        figure = pyplot.figure(
            figsize=(
                row_label_width + cell_width * num_cols,
                col_label_height + cell_height * current_figure_rows))

        label_height_ratio = col_label_height / float(cell_height)
        label_width_ratio = row_label_width / float(cell_width)
        gridspec = pyplot.GridSpec(
            current_figure_rows + 1,
            num_cols + 1,
            height_ratios=([label_height_ratio] +
                [1] * current_figure_rows),
            width_ratios=[label_width_ratio] + [1] * num_cols)

        # The first row of each figure gives the column labels.
        for col_num in range(num_cols):
            ax = pyplot.subplot(gridspec.new_subplotspec((0, col_num + 1)))
            ax.axis("off")
            ax.text(**text_args(
                col_labels[col_num],
                x=0.5,
                y=1.0,
                rotation=90,
                fontsize='x-large'))
            
        # Subsequent rows are the variants.
        for row_num in range(current_figure_row_start, current_figure_row_end):
            local_row_num = row_num - current_figure_row_start

            # First column gives row label.
            ax_row_label = pyplot.subplot(
                gridspec.new_subplotspec((local_row_num + 1, 0)))
            ax_row_label.text(**text_args(
                row_labels[row_num],
                x=0.5,
                y=0.5,
                transform=ax_row_label.transAxes,
                verticalalignment="center"))
            ax_row_label.axis("off")

            # List of (display name, [col1 count, col2 count, ...]) pairs.
            allele_vectors = allele_vectors_matrix[row_num]
            display_alleles = [allele for (allele, _) in allele_vectors]
            col_to_counts_vector = numpy.array(
                [vector for (_, vector) in allele_vectors]).T
            axis_with_a_pie_chart = None
            for col_num in range(num_cols):
                ax = pyplot.subplot(
                    gridspec.new_subplotspec((local_row_num + 1, col_num + 1)))
                ax.axis("off")

                vector = col_to_counts_vector[col_num]

                radius = radii[row_num][col_num]
                if radius > 0:
                    ax.pie(vector, radius=radius, colors=pie_chart_colors)
                    axis_with_a_pie_chart = ax

                ax.text(**text_args(
                    cell_labels[row_num][col_num],
                    x=0.5,
                    y=-0.2,
                    transform=ax.transAxes,
                    fontsize="large",
                    horizontalalignment="center"))

            if axis_with_a_pie_chart:
                axis_with_a_pie_chart.legend(
                    display_alleles,
                    bbox_transform=ax_row_label.transAxes,
                    **legend_properties)

        pyplot.subplots_adjust(**subplots_adjust_properties)
        yield figure
