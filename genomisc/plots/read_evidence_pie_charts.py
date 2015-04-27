import math
from collections import Counter
import logging

from matplotlib import pyplot

VARIANT_LABEL = "\n".join([
    "#{variant_num}",
    "{variant.locus.contig}:{variant.locus.position}",
    "{gene_names}",
    "{variant.ref}->{variant.alt}",
    "{effect}",
    "{extra}"
])

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

def plot(
        variants,
        sample_names,
        evidence_function,
        min_percent_to_show_allele=1.0,
        max_alleles=6,
        variants_per_figure=10,
        sample_name_label_height=5.0,
        variant_label_width=3.0,
        cell_width=2.0,
        cell_height=2.5,
        sample_display_name_function=lambda name: name,
        sample_name_text_extra_properties={},
        variant_label_string=VARIANT_LABEL,
        variant_label_extra_function=lambda variant: "",
        variant_label_text_extra_properties={},
        pie_chart_colors=PIE_CHART_COLORS,
        cell_label_extra_function=lambda variant, sample: "",
        cell_label_text_extra_properties_function=lambda variant, sample: {},
        radius_function=lambda total: math.pow(total / 1300., 1 / 3.0),
        abbreviate_allele_function=abbreviate_allele,
        legend_extra_properties={},
        subplots_adjust_properties={'hspace': 0.1}):

    sample_display_names = dict(
        (name, sample_display_name_function(name))
        for name in sample_names)

    sample_name_text_properties = {
        'rotation': 90,
        'fontsize': 'x-large',
    }
    sample_name_text_properties.update(sample_name_text_extra_properties)

    variant_label_text_properties = {
        'verticalalignment': "center",
    }
    variant_label_text_properties.update(variant_label_text_extra_properties)

    legend_properties = {
        'bbox_to_anchor': (.15, .8),
        'fontsize': 'large',
        'loc': 'upper left',
    }
    legend_properties.update(legend_extra_properties)
    
    variants_queue = list(variants)
    figures = []
    while variants_queue:
        logging.info("Creating figure %d" % (len(figures) + 1))
        plot_variants = [
            variants_queue.pop(0)
            for _
            in range(min(variants_per_figure, len(variants_queue)))
        ]
        figure = pyplot.figure(
            figsize=(
                variant_label_width + cell_width * len(sample_names),
                sample_name_label_height + cell_height * len(plot_variants)))
        figures.append(figure)

        label_height_ratio = sample_name_label_height / float(cell_height)
        label_width_ratio = variant_label_width / float(cell_width)
        gridspec = pyplot.GridSpec(
            len(plot_variants) + 1,
            len(sample_names) + 1,
            height_ratios=[label_height_ratio] + [1] * len(plot_variants),
            width_ratios=[label_width_ratio] + [1] * len(sample_names))

        # The first row of each figure gives the sample names.
        for (sample_i, sample) in enumerate(sample_names):
            ax = pyplot.subplot(gridspec.new_subplotspec((0, sample_i + 1)))
            ax.axis("off")
            ax.text(
                .5,
                1,
                sample_display_names[sample],
                **sample_name_text_properties)

        # Subsequent rows are the variants.
        for (variant_num, variant) in enumerate(plot_variants):
            # First column gives info about the variant.
            variant_label = variant_label_string.format(
                variant_num=variant_num,
                variant=variant,
                gene_names=" ".join(variant.varcode.gene_names()),
                effect=variant.varcode.top_effect().short_description(),
                extra=variant_label_extra_function(variant)).strip()
            
            ax_variant_labels = pyplot.subplot(
                gridspec.new_subplotspec((variant_num + 1, 0)))
            ax_variant_labels.text(
                0.5,
                .5,
                variant_label,
                transform=ax_variant_labels.transAxes,
                **variant_label_text_properties)
            ax_variant_labels.axis("off")

            # Collect the evidence.
            # "Evidence" means an allele -> count Counter.
            sample_to_evidence = dict(
                (sample, Counter(evidence_function(sample, variant)))
                for sample
                in sample_names)

            total_evidence = Counter()
            alleles_with_sufficient_evidence = set()
            for evidence in sample_to_evidence.values():
                total_evidence += evidence
                total = sum(count for (_, count) in evidence.most_common())
                alleles_with_sufficient_evidence.update(
                    allele for (allele, count)
                    in evidence.most_common()
                    if count * 100.0 / total >= min_percent_to_show_allele)

            display_alleles = [variant.ref, variant.alt] + [
                allele for (allele, _) in total_evidence.most_common()
                if allele not in (variant.ref, variant.alt) and (
                    allele in alleles_with_sufficient_evidence)
            ][:max_alleles]

            undisplayed_alleles = set(
                allele for allele in total_evidence
                if allele not in display_alleles)

            if undisplayed_alleles:
                undisplayed_alleles.add(display_alleles.pop(-1))
                display_alleles.append(None)
                for e in [total_evidence] + list(sample_to_evidence.values()):
                    e[None] = sum(e[allele] for allele in undisplayed_alleles)

            axis_with_a_pie_chart = None
            for (sample_i, sample) in enumerate(sample_names):
                ax = pyplot.subplot(
                    gridspec.new_subplotspec((variant_num + 1, sample_i + 1)))
                ax.axis("off")

                values = [
                    sample_to_evidence[sample][allele]
                    for allele in display_alleles
                ]
                values_label = "\n".join([
                    " ".join(str(d) for d in values),
                    "of %d" % sum(values),
                    cell_label_extra_function(variant, sample)]).strip()

                cell_label_font_properties = {
                    'fontsize': 'large',
                    'horizontalalignment': 'center',
                }
                cell_label_font_properties.update(
                    cell_label_text_extra_properties_function(variant, sample))
                ax.text(
                    .5,
                    0,
                    values_label,
                    transform=ax.transAxes,
                    **cell_label_font_properties)

                radius = radius_function(sum(values))
                if radius > 0:
                    ax.pie(values, radius=radius, colors=pie_chart_colors)
                    axis_with_a_pie_chart = ax

            if axis_with_a_pie_chart:
                axis_with_a_pie_chart.legend(
                    [abbreviate_allele(allele) for allele in display_alleles],
                    bbox_transform=ax_variant_labels.transAxes,
                    **legend_properties)

        pyplot.subplots_adjust(**subplots_adjust_properties)

    return figures
