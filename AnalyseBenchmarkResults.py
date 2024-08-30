"""Analysis script dedicated to the HADDOCK3 + haddock-runner benchmarks.

Generates multiple-plots analysing different scenarios performances.
- Barplots: Standard best performing model among top X from all targets.
- Melquiplots: Per-target complex performances among top 200.
- Violinplots: Performance distribution among top X from all targets.

Please modify the Global variable: CAPRIEVAL_STEPS to suite your needs.
It is used to generate nice title to the caprieval steps.

Usage:
>python3 AnalyseBenchmarkResults.py <path/to/benchmark/dir/to/analyse/>
"""

import argparse
import glob
import json
import os
import sys
import zipfile

from pathlib import Path
from typing import Callable, Optional, Union

# Try to load external libraries
try:
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
except ModuleNotFoundError:
    sys.exit(
        "\n[ERROR]: Issue when loading the external libraries.\n\n"
        "This script is using `numpy` and `matplotlib` libraries.\n"
        "Please make sure they are accessible with current environment.\n"
        "(e.g.: >pip install numpy matplotlib)\n"
        )


__version__ = "1.0.5"  # August 2024
__author__ = ", ".join((
    "BonvinLab",
    "Computational Structural Biology group",
    "Utrecht University",
    "the Netherlands",
    "Europe",
    "Planet Earth",
    "Milky Way",
    ))
__dev__ = (
    "Victor G.P. Reys",
    )


####################
# GLOBAL VARIABLES #
####################
# Define custom caprieval steps names
# NOTE: "Feel free to modify this `CAPRIEVAL_STEPS` dict content to fit your experiment"  # noqa : E501
# This dict must have:
#  as keys   -> Index of the caprieval stage (used to parse data)
#  as values -> Name to give to this stage (used as legends plots)
# NOTE: e.g.: for the following run: [topoaa, rigidbody, caprieval, seletop, caprieval.1, flexref, caprieval.2, emref, caprieval.3, clustfcc, seletopclusts, caprieval.4]  # noqa : E501
# CAPRIEVAL_STEPS = {
    # '02': 'rigidbody',
    # '04': 'seletop 200',
    # '06': 'flexref',
    # '08': 'emref',
    # '11': 'top 4 models per fcc clust',
    # }
CAPRIEVAL_STEPS = {
    '02': 'rigidbody',
    '04': 'seletop 200',
    '06': 'flexref',
    '08': 'emref',
    '11': 'top 4 models per fcc clust',
    }

# Set threshold of top X structures to take into account
TOP_X_THRESHOLDS = (1, 5, 10, 20, 50, 100, 200, 500, 1000, )
# Set number of entries to display in melquiplot
MELQUIPLOT_NB_ENTRIES = 200

# CAPRI performances classes
# NOTE: for each class, we define the lower and upper limit
ALL_PERFORMANCES_CLASSES = {
    "protein": {
        "irmsd": {
            "High": (0, 1),
            "Medium": (1, 2),
            "Acceptable": (2, 4),
            "Near-acceptable": (4, 6),
            "Low": (6, 99999),
            "Missing": (-2, -0.5),
            },
        "dockq": {
            "High": (0.8, 1),
            "Medium": (0.6, 0.8),
            "Acceptable": (0.5, 0.6),
            "Near-acceptable": (0.4, 0.5),
            "Low": (0, 0.4),
            "Missing": (-2, -0.5),
            },
        },
    "peptide": {
        "irmsd": {
            "High": (0, 0.5),
            "Medium": (0.5, 1),
            "Acceptable": (1, 2),
            "Near-acceptable": (2, 3),
            "Low": (3, 99999),
            "Missing": (-2, -0.5),
            },
        "dockq": {
            "High": (0.895, 1),
            "Medium": (0.71, 0.895),
            "Acceptable": (0.43, 0.71),
            "Near-acceptable": (0.35, 0.43),
            "Low": (0, 0.35),
            "Missing": (-2, -0.5),
            },
        },
    # FIXME: Optimize irmsd and dockq values for glycan
    "glycan": {
        "irmsd": {
            "High": (0, 0.5),
            "Medium": (0.5, 1),
            "Acceptable": (1, 2),
            "Near-acceptable": (2, 3),
            "Low": (3, 99999),
            "Missing": (-2, -0.5),
            },
        "dockq": {
            "High": (0.895, 1),
            "Medium": (0.71, 0.895),
            "Acceptable": (0.43, 0.71),
            "Near-acceptable": (0.35, 0.43),
            "Low": (0, 0.35),
            "Missing": (-2, -0.5),
            },
        "ilrmsd": {
            "High": (0, 1),
            "Medium": (1, 2),
            "Acceptable": (2, 3),
            "Near-acceptable": (3, 4),
            "Low": (4, 99999),
            "Missing": (-2, -0.5),
            },
        },
    }
PERFORMANCES_CLASSES = ALL_PERFORMANCES_CLASSES["protein"]

# Add color mapper
COLORS_MAPPER = {
    "High": "darkgreen",
    "Medium": "lightgreen",
    "Acceptable": "lightblue",
    "Near-acceptable": "gainsboro",
    "Low": "white",
    "Missing": "dimgrey",
    }

# Performance order
PERF_ORDER = (
    "High",
    "Medium",
    "Acceptable",
    "Near-acceptable",
    "Low",
    "Missing",
    )
# DPI of the generated figures
DPI = 400


####################
# DEFINE FUNCTIONS #
####################
def gen_graph(
        ax: plt.Axes,
        high: list,
        med: list,
        acc: list,
        nacc: list,
        low: list,
        miss: list,
        top: list,
        width: float = 0.5,
        percentage: bool = True,
        ) -> None:
    """Plot a barplot on the provided axis `ax`.

    Parameters
    ----------
    ax : `matplotlib.pyplot.Axes`
        The axis on which to draw the plot.
    high : list
        List containing number of `high performances` models
        for each threshold in `top`.
    med : list
        List containing number of `medium performances` models
        for each threshold in `top`.
    acc : list
        List containing number of `acceptable performances` models
        for each threshold in `top`.
    nacc : list
        List containing number of `near acceptable performances` models
        for each threshold in `top`.
    low : list
        List containing number of `low performances` models
        for each threshold in `top`.
    miss : list
        List containing number of missing model data
        for each threshold in `top`.
    top : list
        List of number of entries take into consideration.
    width : float
        Width of the bar to draw.
    percentage : bool
        If true, number of entries are converted into % sucess
    """
    # Performances mapper
    performances = {
        "High": high,
        "Medium": med,
        "Acceptable": acc,
        "Near-acceptable": nacc,
        "Low": low,
        "Missing": miss,
        }
    if percentage:
        # Compute total at each position
        indices_total: dict[int, int] = {}
        for label, counts in performances.items():
            for i, val in enumerate(counts):
                total = indices_total.setdefault(i, 0)
                indices_total[i] = total + val
        # Compute percentages
        percentage_perfs: dict[str, list[float]] = {}
        for label, counts in performances.items():
            for ind, val in enumerate(counts):
                try:
                    percent = 100 * val / indices_total[ind]
                except ZeroDivisionError:
                    percent = 0
                finally:
                    precent_list = percentage_perfs.setdefault(label, [])
                    precent_list.append(percent)
        performances = percentage_perfs
    # X labels
    tops = [f'Top{v}' for v in top]
    # initialize first bottom values with 0s
    bottom = np.zeros(len(high))

    # Loop over performances
    for label in PERF_ORDER:
        # point performances data
        perfs = performances[label]
        # draw bars
        ax.bar(
            tops,
            perfs,
            width,
            label=label,
            bottom=bottom,
            color=COLORS_MAPPER[label],
            )
        # increment bottom value
        bottom += perfs
    # New labels for percentage specific displaying
    if percentage:
        yticks = [0, 25, 50, 75, 100]
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)
        # Draw horizontal lines for better reading
        xlims = ax.get_xlim()
        xstart = np.floor(xlims[0])
        xend = np.ceil(xlims[1])
        for yposition in yticks[1:-1]:
            ax.plot(
                [xstart, xend],
                [yposition, yposition],
                linestyle="dashed",
                color="gray",
                alpha=0.7,
                )
    # Orient X labels
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right')
    ylabel = "Nb. entries" if not percentage else "% sucess rate"
    ax.set_ylabel(ylabel)


def clear_plt() -> None:
    """Clear all previous instances/data generated by matplotlib."""
    plt.gca()
    plt.cla()
    plt.clf()


def gen_violin(
        ax: plt.Axes,
        perf_data: list,
        labels: list,
        metric: str = "",
        ) -> None:
    """Draw a violinplot on the axis.

    inspired from:
    https://matplotlib.org/stable/gallery/statistics/customized_violin.html

    Parameters
    ----------
    ax : `matplotlib.pyplot.Axes`
        The axis on which to draw the plot.
    perf_data : list
        List of performances values.
    labels : list
        List of labels (same order as `perf_data`).
    """
    # Draw it
    parts = ax.violinplot(
        perf_data,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        )

    # Get color ramp
    nbcolors = 10 if len(perf_data) <= 10 else 20
    colorramp = mpl.colormaps[f'tab{nbcolors}']

    # Modify colors
    for vi, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colorramp((vi + 0.5) / nbcolors))
        pc.set_edgecolor('black')
        pc.set_alpha(1)
    quartile1, medians, quartile3 = np.percentile(
        perf_data,
        [25, 50, 75],
        axis=1,
        )
    whiskers = np.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(perf_data, quartile1, quartile3)
        ])
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
    inds = np.arange(1, len(medians) + 1)
    ax.scatter(inds, medians, marker='_', color='white', s=30, zorder=3)
    ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
    ax.vlines(inds, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)
    ax.set_ylabel(metric)

    # Set labels
    if labels:
        ax.set_xticks(list(range(1, len(labels) + 1)))
        ax.set_xticklabels(
            labels,
            rotation=40,
            ha='right',
            rotation_mode='anchor',
            )


def adjacent_values(
        vals: list[float],
        q1: float,
        q3: float,
        ) -> tuple[float, float]:
    """Find adjacent values.

    Inspired from:
    https://matplotlib.org/stable/gallery/statistics/customized_violin.html

    Parameters
    ----------
    vals : list
        List of values
    q1 : float
        Value of the first quartil
    q3 : float
        Value of the 3rd quartil

    Return
    ------
    lower_adjacent_value : float
        Closest true value under q1
    upper_adjacent_value : float
        Closest true value above q3
    """
    # Finds closest true value above q3
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])
    # Find closest true value under q1
    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def stage_name(cname: str) -> str:
    """Try to return the user defined stage name, or return default.

    Parameters
    ----------
    cname : str
        Index of the caprieval stage

    Returns
    -------
    name : str
        Name of the stage
    """
    try:
        name = CAPRIEVAL_STEPS[cname]
    except KeyError:
        name = f"{cname}_caprieval"
    return name


def gen_full_comparison_violins(
        scenars_perfs: dict,
        basepath: str = "./",
        title: str = "",
        metric: str = "",
        progress: bool = True,
        ) -> None:
    """Combine all scenarios caprieval within same plot.

    Parameters
    ----------
    scenars_perfs : dict
        Dictionary of all scenario stages performances.
    basepath : str
        Basepath where to write plot.
    title : str
        Title of the figure.
    """
    # Clear pervious instances of matplotlib
    clear_plt()
    # Compute number of rows (scenarios)
    scenars_order = sorted(scenars_perfs)
    nb_scenar = len(scenars_order)
    # Compute number of colums (capri steps)
    steps_order = sorted(scenars_perfs[scenars_order[0]])
    nb_steps = len(steps_order)
    # Get number of threshodls
    tops_order = sorted(
        scenars_perfs[scenars_order[0]][steps_order[0]]['values'],
        )
    nb_thresh = len(tops_order)
     # Compute total number of plots
    total_plots = nb_thresh * nb_steps
    processed = 0

    # Initate figures / axis
    fig, axes = plt.subplots(
        figsize=((nb_steps * 4) + 1, (nb_thresh * 3) + 1),
        nrows=nb_thresh,
        ncols=nb_steps,
        sharey=True,
        sharex=True,
        )

    # Loop over rows
    for ri, topx in enumerate(tops_order):
        for ci, cname in enumerate(steps_order):
            processed += 1
            if progress:
                print(f"{100 * processed / total_plots:>6.2f} %", end="\r")
            # Point axis
            ax = axes[ri][ci]
            # Build sub-dt dict
            perf_data = [
                sorted(scenars_perfs[scenar][cname]['values'][topx])
                for scenar in scenars_order
                ]
            # Set labels on last row only
            labels = None
            if ri + 1 == nb_thresh:
                labels = [so.replace('scenario-', '') for so in scenars_order]
            # Write bars
            gen_violin(
                ax,
                perf_data,
                labels,
                metric=metric,
                )

    # Add columns titles
    pad = 5
    for ax, cname in zip(axes[0], steps_order):
        # Annotate column
        ax.annotate(
            stage_name(cname),
            xy=(0.5, 1),
            xytext=(0, pad),
            xycoords='axes fraction',
            textcoords='offset points',
            size='large',
            ha='center',
            va='baseline',
            )
    # Add rows titles
    for ax, topx in zip(axes[:, 0], tops_order):
        ax.annotate(
            f'Top {topx}',
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad - pad, 0),
            xycoords=ax.yaxis.label,
            textcoords='offset points',
            size='large',
            ha='right',
            va='center',
            )

    # Add figure title
    fig.suptitle(title, fontsize=16)

    # Get color ramp
    nbcolors = 10 if nb_scenar <= 10 else 20
    colorramp = mpl.colormaps[f'tab{nbcolors}']
    # Add bars legend
    fig.legend(
        [
            mpatches.FancyBboxPatch(
                (-0.025, -0.05), 0.05, 0.1, ec="none",
                boxstyle=mpatches.BoxStyle("Round", pad=0.02),
                color=colorramp((si + 0.5) / nbcolors),
                )
            for si, _perfclass in enumerate(scenars_order)
            ],
        [so.replace('scenario-', '') for so in scenars_order],
        loc='outside lower center',
        ncols=4,
        title="Screnarios",
        )

    plt.gca().set_ylim(bottom=0)

    # adjust border to let annotations fit inside graph
    fig.subplots_adjust(left=0.08, top=0.95, bottom=0.12, right=0.98)

    # save figure
    plt.savefig(f"{basepath}_violins.png", format='png', dpi=DPI)
    return


def gen_full_comparison_melquiplots(
        scenars_perfs: dict,
        perf_dtype: str = "irmsd",
        basepath: str = "./",
        progress: bool = True,
        ) -> str:
    """Generate multiple melquiplots for each scenario.

    Parameters
    ----------
    scenars_perfs : dict
        Dict containing performances for each scenario.
    perf_dtype : str, optional
        Model quality metric to use, by default "irmsd"
    title : str, optional
        Prefix to give to archive, by default "benchmark_melquis"
    basepath : str, optional
        Where to write the files, by default "./"

    Returns
    -------
    archive_path : str
        Path to the generated archive.zip
    """
    # Clear pervious instances of matplotlib
    clear_plt()
    
    # Set progression variables
    all_generated_melquis: list[str] = []
    processed = 0
    total_plots = len(scenars_perfs)

    # Loop over scenarios
    for scenar_name, scenar_perfs in scenars_perfs.items():
        processed += 1
        if progress:
            print(f"{100 * processed / total_plots:>6.2f} %", end="\r")
        # Generate melquiplot for this scenario
        scenar_melqui_path = make_scenar_melquiplots(
            scenar_perfs,
            perf_dtype=perf_dtype,
            title=scenar_name,
            basepath=basepath,
            )
        all_generated_melquis.append(scenar_melqui_path)

    # Generate archive of melqui plots
    archive_path = Path(f"{basepath}_melquiplots.zip")
    path = archive_path.parent
    archive_name = archive_path.name
    initdir = os.getcwd()
    os.chdir(path)
    with zipfile.ZipFile(archive_name, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for figure in all_generated_melquis:
            zipf.write(Path(figure).name)
    os.chdir(initdir)

    # Remove all original files
    for generated_melqui in all_generated_melquis:
        os.remove(generated_melqui)

    # Return archive path
    return archive_path


def make_scenar_melquiplots(
        scenar_perfs: dict,
        perf_dtype: str = "irmsd",
        title: str = "melquiplot",
        basepath: str = "./",
        ) -> str:
    """Generate multiples melquiplots for each Caprieval steps of a scenario.

    Parameters
    ----------
    scenar_perfs : dict
        Dict containing performances of a scenario.
    perf_dtype : str, optional
        Model quality metric to use, by default "irmsd"
    title : str, optional
        Title to the figure, by default "melquiplot"
    basepath : str, optional
        Where to write the files, by default "./"

    Returns
    -------
    figpath : str
        Path to the generated figure.
    """
    # Count nb_steps and entries
    steps = sorted(scenar_perfs)
    nb_rows = len(steps)
    nb_entries = len(scenar_perfs[steps[0]].keys())
    dtype_perf_classes = PERFORMANCES_CLASSES[perf_dtype]
    # Initate figures / axis
    fig, axes = plt.subplots(
        figsize=((nb_entries * 1) + 1, (nb_rows * 5) + 3),
        nrows=nb_rows,
        ncols=1,
        )
    # Loop over stages
    for si, (stage, stage_perfs) in enumerate(scenar_perfs.items()):
        # Point axis
        ax = axes[si]
        # Draw a melquiplot on this axis
        melquiplot(ax, stage_perfs, width=min(1, 5 / nb_entries))
        # Add title to graph
        ax.annotate(
            stage_name(stage),
            xy=(0.5, 1),
            xytext=(0, 5),
            xycoords='axes fraction',
            textcoords='offset points',
            size='large',
            ha='center',
            va='baseline',
            )
    # Add figure title
    fig.suptitle(title, fontsize=16)
    # Add Legend
    legend_data = [
        (
            plt.Rectangle((0, 0), 1, 1, fc=COLORS_MAPPER[perfclass]),
            rf'{perfclass} | {dtype_perf_classes[perfclass][0]} <= {perf_dtype.upper()} < {dtype_perf_classes[perfclass][1]}$\AA$',  # noqa : E501
            )
        for perfclass in PERF_ORDER
        ]
    legend_proxies, legend_labels = zip(*legend_data)
    # Add bars legend
    fig.legend(
        legend_proxies,
        legend_labels,
        loc='outside lower center',
        ncols=len(PERF_ORDER),
        title="performance classes",
        )

    # adjust border to let annotations fit inside graph
    fig.subplots_adjust(left=0.02, top=0.95, bottom=0.06, right=0.98)

    # save figure
    figpath = f"{basepath}_{title}_melquiplot.png"
    plt.savefig(figpath, format='png', dpi=DPI)
    return figpath


def melquiplot(
        ax: plt.Axes,
        pdb_perfs: dict[str, list[str]],
        width: float = 0.3,
        ) -> None:
    """Draw a melquiplot on an sub-figure with provided input data.

    Parameters
    ----------
    ax : plt.Axes
        The axis on which to draw the Melquiplot
    pdb_perfs : dict[str, list[str]]
        Performances for each entry at a give stage.
    """
    width = 0.1
    max_stack_y = 0
    pdb_labels = sorted(pdb_perfs)
    stack_h_labels_pos: list[float] = []
    # Loop over each input entry
    for entry_index, pdb in enumerate(pdb_labels, start=0):
        # Point data
        perfs = pdb_perfs[pdb][MELQUIPLOT_NB_ENTRIES]
        x_coord = entry_index * width
        y_coord = 0
        # Loop over perfs
        for perf_label in perfs:
            # Draw a bar
            ax.bar(
                x_coord,
                2,
                align="edge",
                bottom=y_coord,
                color=COLORS_MAPPER[perf_label],
                width=width * 0.98,
                edgecolor='none',
                linewidth=0,
                )
            y_coord += 1

        max_stack_y = max(max_stack_y, y_coord)
        stack_h_labels_pos.append(x_coord + (width / 2))

    # Aesthetics
    ax.set_xlim((0, len(pdb_labels) * width))
    ax.set_ylim((0, max_stack_y))
    ax.set_ylabel('Energy Score Ranking (lower is better)')
    ax.set_xticks(stack_h_labels_pos)
    ax.set_xticklabels(pdb_labels, rotation=45, fontsize='small')
    ax.xaxis.set_ticks_position('none')


def melquiplot_original(
        ax: plt.Axes,
        pdb_perfs: dict[str, list[str]],
        ) -> None:
    """Draw a melquiplot on an sub-figure with provided input data.

    Parameters
    ----------
    ax : plt.Axes
        The axis on which to draw the Melquiplot
    pdb_perfs : dict[str, list[str]]
        Performances for each entry at a give stage.
    """
    max_stack_v = 0
    pdb_labels = sorted(pdb_perfs)
    stack_h_labels_pos: list[float] = []
    # Loop over each input entry
    for entry_index, pdb in enumerate(pdb_labels, start=1):
        # Point data
        perfs = pdb_perfs[pdb][MELQUIPLOT_NB_ENTRIES]
        stack_v = 0
        # Loop over perfs
        for perf_label in perfs:
            # Draw a bar
            ax.bar(
                entry_index - 0.5,
                2,
                bottom=stack_v,
                color=COLORS_MAPPER[perf_label],
                width=0.99,
                edgecolor='none',
                linewidth=0,
                )
            stack_v += 1

        max_stack_v = max(max_stack_v, stack_v)
        stack_h_labels_pos.append(entry_index - 0.501)

    # Aesthetics
    ax.set_xlim((0.05, len(pdb_labels) + 0.05))
    ax.set_ylim((0, max_stack_v))
    ax.set_ylabel('Energy Score Ranking (lower is better)')
    ax.set_xticks(stack_h_labels_pos)
    ax.set_xticklabels(pdb_labels, rotation=45, fontsize='small')
    ax.xaxis.set_ticks_position('none')


def gen_full_comparison_barplots(
        scenars_perfs: dict,
        basepath: str = "./",
        title: str = "",
        progress: bool = True,
        no_percentage: bool = False,
        ) -> None:
    """Combine all scenarios caprieval within same plot.

    Parameters
    ----------
    scenars_perfs : dict
        Dictionary of all scenario stages performances.
    basepath : str
        Basepath where to write plot.
    title : str
        Title of the figure.
    """
    # Clear pervious instances of matplotlib
    clear_plt()
    # Compute number of rows
    rows_order = sorted(scenars_perfs)
    nb_rows = len(rows_order)
    # Compute number of colums
    cols_order = sorted(scenars_perfs[rows_order[0]])
    nb_cols = len(cols_order)
    # Compute total number of plots
    total_plots = nb_rows * nb_cols
    processed = 0

    # Initate figures / axis
    fig, axes = plt.subplots(
        figsize=((nb_cols * 4) + 1, (nb_rows * 3) + 4),
        nrows=nb_rows, ncols=nb_cols,
        )

    # Loop over rows
    for ri, rname in enumerate(rows_order):
        # Point row axe(s)
        if len(rows_order) == 1:
            axr = axes
        else:
            axr = axes[ri]
        for ci, cname in enumerate(cols_order):
            # Display progression
            processed += 1
            if progress:
                print(f"{100 * processed / total_plots:>6.2f} %", end="\r")
            # Point column axe(s)
            if len(cols_order) == 1:
                ax = axr
            else:
                ax = axr[ci]
            # Point data
            perf_data = scenars_perfs[rname][cname]['classes']
            # Get sorted entreis
            topx = sorted(perf_data)
            # attribute perf classes at each topX model
            high = [perf_data[top]['High'] for top in topx]
            med = [perf_data[top]['Medium'] for top in topx]
            acc = [perf_data[top]['Acceptable'] for top in topx]
            nacc = [perf_data[top]['Near-acceptable'] for top in topx]
            low = [perf_data[top]['Low'] for top in topx]
            miss = [perf_data[top]['Missing'] for top in topx]
            # Write bars
            gen_graph(
                ax,
                high,
                med,
                acc,
                nacc,
                low,
                miss,
                topx,
                percentage=not no_percentage,
                )

    # Set padding between two plots
    pad = 5
    # Search for first row
    if len(rows_order) == 1:
        first_row = axes
    else:
        first_row = axes[0]

    # Add columns titles
    for ax, cname in zip(first_row, cols_order):
        # Add column name
        ax.annotate(
            stage_name(cname),
            xy=(0.5, 1),
            xytext=(0, pad),
            xycoords='axes fraction',
            textcoords='offset points',
            size='large',
            ha='center',
            va='baseline',
            )

    # Find all first row columns
    if len(rows_order) == 1:
        axrows = [axes]
    else:
        axrows = axes
    if len(cols_order) == 1:
        axrows_firstcols = axrows
    else:
        axrows_firstcols = [ax[0] for ax in axrows]

    # Add rows titles
    for ax, rname in zip(axrows_firstcols, rows_order):
        ax.annotate(
            rname.replace('scenario-', ''),
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad - pad, 0),
            xycoords=ax.yaxis.label,
            textcoords='offset points',
            size='large',
            ha='right',
            va='center',
            )

    # Add figure title
    fig.suptitle(title, fontsize=16)

    # Add bars legend
    fig.legend(
        [mpatches.FancyBboxPatch(
            (-0.025, -0.05), 0.05, 0.1, ec="none",
            boxstyle=mpatches.BoxStyle("Round", pad=0.02),
            color=COLORS_MAPPER[perfclass],
            )
         for perfclass in PERF_ORDER],
        PERF_ORDER,
        loc='outside lower center',
        ncols=len(PERF_ORDER),
        title="performance classes",
        )

    # Adjust border to let annotations fit inside graph
    fig.subplots_adjust(left=0.15, top=0.9, bottom=0.15, right=0.98)

    # Save figure
    plt.savefig(f"{basepath}_capribarpolots.png", format='png', dpi=DPI)
    return


def get_pdb_entries(basepath: str) -> list:
    """Retrieve list of PDBid.

    Parameters
    ----------
    basepath : str
        Path to a scenario directory containing pdb entries.

    Return
    ------
    pdbids : str
        List of pdb entries benchmarked in this scenario.
    """
    pdbids = [
        pdbid_path.stem
        for pdbid_path in Path(basepath).glob("*/")
        if pdbid_path.is_dir()
        ]
    return pdbids


def get_scenarios_names(basepath: str) -> list:
    """Retrieve list of scenario names.

    Parameters
    ----------
    basepath : str
        Path to the benchmark directory to analyse containing X scenarios.

    Return
    ------
    scenarios_names : list
        List of tested scenario names.
    """
    scenario_paths = glob.glob(f'{basepath}scenario*/')
    scenarios_names = [sp.split('/')[-2] for sp in scenario_paths]
    return scenarios_names


def scenario_name_2_threshold(scenar_name: str) -> float:
    """Gather threshold value from scenario name.

    NOTE: function used for CPORT-ARCTIC3D-BM5 benchmark

    Parameters
    ----------
    scenar_name : str
        Name of a scenario.

    Return
    ------
    threshold : float
        CPORT threshold value used in this scenario.
    """
    # Get last part of scenario name
    str_thresh = scenar_name.split('-')[-1]
    # Replace underscore by a dot
    dot_thresh = str_thresh.replace('_', '.')
    # Cast it to float
    threshold = float(dot_thresh)
    return threshold


def get_caprieval_stages(basepath: str) -> list[str]:
    """Retrieve all caprieval stages inside a haddock3 run.

    Parameters
    ----------
    basepath : str
        Path to a haddock3 run directory `rundir`.

    Return
    ------
    caprieval_stages : list
        List of caprieval module indexes
    """
    caprieval_paths = glob.glob(f'{basepath}*caprieval/')
    caprieval_stages = [
        hd3_module_2_stage(caprip.split('/')[-2])
        for caprip in caprieval_paths
        ]
    return caprieval_stages


def hd3_module_2_stage(indexed_modulename: str) -> str:
    """Split a haddock3 module name and retrieve it id.

    Parameters
    ----------
    indexed_modulename : str
        Directory name of an indexed haddock3 module name.
        

    Return
    ------
    stage : str
        Index of the haddock3 module.
    """
    stage = indexed_modulename.split('_')[0]
    return stage


def map_data(
        basepath: str,
        subset_scenarios: Optional[list[str]] = None,
        ) -> dict[str, dict[str, dict[str, dict[str, str]]]]:
    """Map data in one analysis dict to accessit easily.

    Parameters
    ----------
    basepath : str
        Path where the benchmarking scenarios can be found.
    subset_scenarios : Optional[list[str]]
        List of scenario names on which to perform the analysis.
   
    Return
    ------
    dtmap : dict[str, dict[str, dict[str, dict[str, str]]]]
        Dictionary mapping scenarios/pdbids to caprieval tsv paths.
        Structure:
        {"scenario_id":
            {"caprieval_stage":
                {"PDBid":
                    {"complexe_id": path/to/caprieval_ss.tsv,
                     ...,
                    }
                }
            }
        }
    """
    # Initiate mapper
    dtmap: dict[str, dict[str, dict[str, dict[str, str]]]] = {}
    # Gather all directories
    pdbids = get_pdb_entries(basepath)

    # Initiate scenarios names holder
    all_scenarios: list[str] = []
    if subset_scenarios:
        all_scenarios += subset_scenarios
    # Search for all available scenarios
    else:
        # Gather all scenarios
        for pdbid in pdbids:
            scenarios = get_scenarios_names(f'{basepath}{pdbid}/')
            all_scenarios += scenarios
        all_scenarios += list(set(all_scenarios))

    # Make sure all pdbs have all scenarios
    for pdbid in pdbids:
        for scenario in all_scenarios:
            scenar_rundir = f"{basepath}{pdbid}/{scenario}/run1/"
            assert os.path.exists(scenar_rundir), \
                f"[ERROR] could not find scenario `{scenario}` directory for entry `{pdbid}` at: {scenar_rundir}"
   
    all_caprieval_stages = []
    # Add scenario data to data maper
    for scenario in all_scenarios:
        # Loop over pdb ids
        for pdbid in pdbids:
            # Generate scenario basepath
            scenario_bp = f'{basepath}{pdbid}/{scenario}/run1/'
            # Retrieve caprieval stages
            caprieval_stages = get_caprieval_stages(scenario_bp)
            all_caprieval_stages += caprieval_stages
    all_caprieval_stages = sorted(list(set(all_caprieval_stages)))
    
    # Make sure all stages are computed for all pdb in all scenarios...
    for scenario in all_scenarios:
        for pdbid in pdbids:
            for stage in all_caprieval_stages:
                # Build caprieval tsv filepath
                caprieval_tsv_path = f"{basepath}{pdbid}/{scenario}/run1/{stage}_caprieval/capri_ss.tsv"  # noqa : E501
                assert os.path.exists(caprieval_tsv_path), \
                    f"[ERROR] could not access CAPRIEVAL results file at: {caprieval_tsv_path}\n- Stage {stage}\n- Scenario `{scenario}`\n- Target `{pdbid}`"

    # Gather all data
    for scenario in all_scenarios:
        dtmap[scenario] = {}
        for stage in all_caprieval_stages:
            dtmap[scenario][stage] = {}
            for pdbid in pdbids:
                # Build caprieval tsv filepath
                caprieval_tsv_path = f'{basepath}{pdbid}/{scenario}/run1/{stage}_caprieval/capri_ss.tsv'  # noqa : E501
                # Hold datapath
                dtmap[scenario][stage][pdbid] = caprieval_tsv_path

    return dtmap


def analyse_scenario(
        scenario_dt: dict,
        _entries_thresholds: list,
        sort_dtype: str = 'haddock-score',
        perf_dtype: str = 'irmsd',
        ) -> tuple[dict[str, dict], dict[str, dict]]:
    """Process the analysis of a scenario.

    Parameters
    ----------
    scenario_dt : dict
        Dictionary of scenario data containing
        - as keys; index of the caprieval stage
        - as values; Dictionary mapping pdbid to caprieval files
    entries_thresholds : list
        List of entries theshold to take into consideration.
    sort_dtype : str
        Key used to sort data.
    perf_dtype : str
        Key used to define performance.
    output : str
        Basepath where to write the data

    Return
    ------
    scenario_stages_perfs : dict
        Dictionary of stage performances, containing
        - as keys; index of the caprieval stage
        - as values; best performing values and classes at each threshold
    """
    if not MELQUIPLOT_NB_ENTRIES in _entries_thresholds:
        entries_thresholds = [e for e in _entries_thresholds]
        entries_thresholds.append(MELQUIPLOT_NB_ENTRIES)
    else:
        entries_thresholds = _entries_thresholds

    # Initiate all stages performances holder
    scenario_stages_perfs = {}
    scenario_stages_pdb_perfs = {}
    # Loop over caprieval data
    for stage, stage_perfs_mapper in scenario_dt.items():
        # Initiate stage performance holder
        pdb_best_perfs = {}
        pdb_perfs = {}
        # Loop over pdbs
        for pdbid, caprieval_filepath in stage_perfs_mapper.items():
            # Gather performances
            best_perfs_h, top_x_perfs = analyse_caprieval_performances(
                caprieval_filepath,
                entries_thresholds,
                sort_dtype=sort_dtype,
                perf_dtype=perf_dtype,
                )
            # Hold best performances
            pdb_best_perfs[pdbid] = best_perfs_h
            # Hold top X perf classes
            pdb_perfs[pdbid] = {
                topx: [perf_to_class(v, dtype=perf_dtype) for v in perfs]
                for topx, perfs in top_x_perfs.items()
                }

        # Summerize all pdb best performances
        stage_best_perfs = {
            n: [pdb_perf[n] for pdb_perf in pdb_best_perfs.values()]
            for n in entries_thresholds
            }
        # Transform performances values to classes
        stage_class_perfs = {
            n: perfs_to_classes(perfs, dtype=perf_dtype)
            for n, perfs in stage_best_perfs.items()
            }

        # Hold stage performances
        scenario_stages_perfs[stage] = {
            'values': stage_best_perfs,
            'classes': stage_class_perfs,
            }
        scenario_stages_pdb_perfs[stage] = pdb_perfs

    return scenario_stages_perfs, scenario_stages_pdb_perfs


def perfs_to_classes(
        perfs: list,
        dtype: str = "irmsd",
        ) -> dict:
    """Convert performances values into classes.

    Parameters
    ----------
    perfs : list
        List of performance values.
    dtype : str
        Type of the performance.

    Return
    ------
    perfs_classes : dict
        Dictionary holding
        - as keys; names of the classes
        - as values; number of entries falling in this classe
    """
    # Count classes in perfs
    perfs_classes_list = [
        perf_to_class(v, dtype=dtype)
        for v in perfs
        ]
    perfs_classes = {
        cls: perfs_classes_list.count(cls)
        for cls in PERFORMANCES_CLASSES[dtype].keys()
        }
    return perfs_classes


def perf_to_class(value: float, dtype: str = "irmsd") -> str:
    """Cast a performance value into a performance class.

    Parameters
    ----------
    value : float
        Performance value.
    dtype : str
        Name of the data type.

    Return
    ------
    cls : str
        Corresponding class
    """
    # Get dt types classes boundaries
    classes_dt_map = PERFORMANCES_CLASSES[dtype]
    # Get comparison function
    comparison_func = dtype_to_boundary_function(dtype)
    for cls, (lowb, highb) in classes_dt_map.items():
        if comparison_func(lowb, highb, value):
            return cls
    # In case it is not found, return "Missing"
    return PERF_ORDER[-1]


def dtype_to_boundary_function(dtype: str) -> Callable:
    """Get function to make a comparison between two boundaries.

    Parameters
    ----------
    dtype : str
        Name of the data type.

    Return
    ------
    _function : function
        Function used to make a comparison.
    """
    def _include_higher(lowb: float, highb: float, value: float) -> bool:
        return lowb < value <= highb
    def _include_lower(lowb: float, highb: float, value: float) -> bool:
        return lowb <= value < highb
    _function = _include_higher if get_reverse_bool(dtype) else _include_lower
    return _function


def analyse_caprieval_performances(
        capireval_fpath: str,
        entries_thresholds: list,
        sort_dtype: str = 'haddock-score',
        perf_dtype: str = 'irmsd',
        ) -> tuple[dict[int, float], dict[int, list[float]]]:
    """Analyse a CAPRIeval step performances.

    Parameters
    ----------
    capireval_fpath : str
        Path to the caprival file to analyse.
    entries_thresholds : list
        List of entries theshold to take into consideration.
    sort_dtype : str
        Key used to sort data.
    perf_dtype : str
        Key used to define performance.

    Return
    ------
    best_perfs_h : dict
        Dictionary holder best performances at each threshold.
    """
    # Load caprieval data
    caprieval_data = load_caprieval_data(capireval_fpath)
    # Sort entries by perf value
    sorted_complexes = sorted(
        caprieval_data,
        key=lambda k: caprieval_data[k][sort_dtype],
        reverse=get_reverse_bool(sort_dtype),
        )

    # Get performances taking into account increasing number of entries
    best_perfs_h = {
        nb_entries: best_perfs(
            caprieval_data,
            sorted_complexes,
            dtype=perf_dtype,
            n=nb_entries,
            )
        for nb_entries in entries_thresholds
        }
    # Gather top X performances
    top_x_perfs = {
        nb_entries: [
            caprieval_data[entry][perf_dtype]
            for entry in sorted_complexes[:nb_entries]
            ]
        for nb_entries in entries_thresholds
        }

    # Return performances
    return best_perfs_h, top_x_perfs


def best_perfs(
        caprieval_data: dict,
        order: list,
        dtype: str = 'irmsd',
        n: int = 1,
        tolerance: float = 6,
        ) -> float:
    """Obtain best performing value.

    Parameters
    ----------
    caprieval_data : dict
        Dictionary holding performances data for each complexes.
    order : list
        Ordered set of keys in `caprieval_data`.
    dtype : str
        Name of the performance key to take into consideration.
    n : int
        Number of entries to take into consideration.
    tolerance : float
        Percentage of tolerated difference between number of accessible
        entries `len(order)` and number of queried `n` entries.

    Return
    ------
    best_perfomance : float
        Best performing value.
    """
    # Check if length of order << n
    # if (len(order) + (len(order) * tolerance / 100)) < n:
    #     return -1

    # Obtain list of all performances
    all_perfs = [caprieval_data[entry][dtype] for entry in order[:n]]
    # Get best function
    best_function = dtype_to_best_function(dtype)
    # Get best performing value
    best_perfomance = best_function(all_perfs)
    return best_perfomance
        
  
def dtype_to_best_function(dtype: str) -> Callable:
    """Get best function.

    Parameters
    ----------
    dttype : str
        Name of the data type.

    Return
    ------
    _function : function
        Function used to identify best score within list of values.
    """
    _function = max if get_reverse_bool(dtype) else min
    return _function
            

def get_reverse_bool(dt_type: str) -> bool:
    """Check sorting order.

    Parameters
    ----------
    dt_type : str
        Name of the data type.

    Return
    ------
    _reversed : bool
        Weather or not the ordering should be reversed.
    """
    _reversed = True if dt_type in ('fnat', 'dockq', ) else False
    return _reversed


def load_caprieval_data(
        tsvpath: str,
        sep: str = '\t',
        ) -> dict[str, dict[str, float]]:
    """Load caprieval data as dict.

    Parameters
    ----------
    tsvpath : str
        Path to the caprieval_ss.tsv file.
    sep : str
        String used to separate data in the file.

    Return
    ------
    data : dict[str, dict[str, float]]
        Dictionary holding data for each complexes.
    """
    data: dict[str, dict[str, float]] = {}
    with open(tsvpath, 'r') as filin:
        for i, _ in enumerate(filin):
            # Split the line
            s_ = _.strip().split(sep)
            # Gather header names
            if i == 0:
                header = s_
                continue
            # Load data
            complex_name = s_[header.index('model')]
            complex_dt_str = {
                'rank': s_[header.index('caprieval_rank')],
                'haddock-score': s_[header.index('score')],
                'irmsd': s_[header.index('irmsd')],
                'lrmsd': s_[header.index('lrmsd')],
                'fnat': s_[header.index('fnat')],
                'ilrmsd': s_[header.index('ilrmsd')],
                'dockq': s_[header.index('dockq')],
                }
            # Cast all values to float
            complex_dt = {
                k: float(v)
                for k, v in complex_dt_str.items()
                }
            # Hold this guy
            data[complex_name] = complex_dt
    return data


def write_json(data: dict, path: str) -> None:
    """Write a json file containing data.

    Parameters
    ----------
    data : dict / list
        The data to be dumped in the file.
    path : str
        Path to the file to write.
    """
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=4)


def load_json(path: str) -> Union[dict, list]:
    """Load a json file.

    Parameters
    ----------
    path : str
        Path to the .json file to load.

    Return
    ------
    data : dict or list
        Python loaded data within file
    """
    with open(path) as f:
        data = json.load(f)
    return data


def get_data_mapper(
        basepath: str,
        overwrite: bool = False,
        outpath: str = "",
        subset_scenarios: Optional[list[str]] = None,
        ) -> dict:
    """Retrieve or generate data mapper.

    Parameters
    ----------
    basepath : str
        Path to the benchmark directory to analyse
    overwrite : bool
        Weather or not an existing benchmark_mapper be overwritten.
        If False, existing benchmark_mapper will be returned.
    outpath : str
        Base path where to write data
    subset_scenarios : Optional[list[str]]
        List of scenario names on which to perform the analysis.

    Return
    ------
    dtmapper : dict
        Dictionary mapping scenarios/pdbids to caprieval tsv paths.
    """
    # Name of the mapper
    mapper_fpath = f'{outpath}_benchmark_mapper.json'
    # Check for no overwrites
    if not overwrite and os.path.exists(mapper_fpath):
        dtmapper = load_json(mapper_fpath)
        return dtmapper

    # Gather data mapper
    dtmapper = map_data(basepath, subset_scenarios=subset_scenarios)
    # Write it
    write_json(dtmapper, mapper_fpath)
    # Return data
    return dtmapper


def gen_outputdir(outputdir: str) -> None:
    """Make sure output directory is created.

    Parameters
    ----------
    outputdir : str
        Path to directory where to store the results
    """
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)


def set_output_path(
        benchmark_directory: str,
        outputpath: str,
        ) -> tuple[str, str]:
    """Initiate output paths and directory.

    Parameters
    ----------
    benchmark_directory : str
        Path to the benchmark directory to analyse
    outputpath : str
        Path to directory where to store the results

    Return
    ------
    basename : str
        Name of the directory to analyse
    base_outputpath : str
        Base path where to write data
    """
    # Initiate output directory
    gen_outputdir(outputpath)
    # Get basename of analysis
    basename = os.path.basename(os.path.dirname(benchmark_directory))
    # Generate baseoutput path
    base_outputpath = f'{outputpath}{basename}'
    return basename, base_outputpath


def vprint(msg: str, silence: bool) -> None:
    """Print on screen

    Parameters
    ----------
    msg : str
        Message to print
    silence : bool
        If true, do not print
    """
    if not silence:
        print(msg)


#################
# Main function #
#################
def main(
        benchmark_directory: str,
        outputpath: str = 'Analysis/',
        scenarios: Optional[list[str]] = None,
        metric: str = "irmsd",
        silent: bool = False,
        no_percentage: bool = False,
        no_capriplots: bool = False,
        no_violinplots: bool = False,
        no_melquiplots: bool = False,
        ) -> None:
    """Run the analysis procedure.

    Parameters
    ----------
    benchmark_directory : str
        Path to the benchmark directory to analyse
    outputpath : str
        Path to directory where to store the results
    scenarios : Optional[list[str]]
        List of scenario names on which to perform the analysis.
    """
    # Initiate output paths
    basename, base_outputpath = set_output_path(
        benchmark_directory,
        outputpath,
        )
    vprint(f"Setting the output directory path: `{outputpath}`", silent)

    # Gather data mapper
    vprint(f"- Searching for data in `{benchmark_directory}`", silent)
    dtmapper = get_data_mapper(
        benchmark_directory,
        overwrite=True,
        outpath=base_outputpath,
        subset_scenarios=scenarios,
        )

    # Initiate all scenarios perforamces mapper
    vprint(
        f"- Loading data from `{benchmark_directory}` "
        f"for {len(dtmapper)} scenario(s): {', '.join(dtmapper)}",
        silent,
        )
    all_scenar_perfs: dict[str, dict] = {}
    all_scenar_melquis: dict[str, dict] = {}
    # Loop over scenario
    for scenar_name, scenario_dt in dtmapper.items():
        # Analyse this scenario
        scenar_best_perfs, scenar_pdb_perfs = analyse_scenario(
            scenario_dt,
            TOP_X_THRESHOLDS,
            sort_dtype='haddock-score',
            perf_dtype=metric,
            )
        # Hold perforamnces
        all_scenar_perfs[scenar_name] = scenar_best_perfs
        all_scenar_melquis[scenar_name] = scenar_pdb_perfs

    # Write data as json
    write_json(all_scenar_perfs, f'{base_outputpath}_performances.json')

    # Draw general graph
    if not no_capriplots:
        vprint("- Generating Bar plots", silent)
        gen_full_comparison_barplots(
            all_scenar_perfs,
            basepath=base_outputpath,
            title=basename,
            progress=not silent,
            no_percentage=no_percentage,
            )

    if not no_violinplots:
        vprint("- Generating Violin plots", silent)
        gen_full_comparison_violins(
            all_scenar_perfs,
            basepath=base_outputpath,
            title=basename,
            metric=metric,
            progress=not silent,
            )

    if not no_melquiplots:
       vprint("- Generating Melqui plots", silent)
       gen_full_comparison_melquiplots(
            all_scenar_melquis,
            perf_dtype=metric,
            basepath=base_outputpath,
            progress=not silent,
            )


###################################
# COMMAND LINE ARGUMENTS HANDLERS #
###################################
def must_end_with_slash(dirpath: str) -> str:
    """Make sure a directory paths is terminating by '/'.
    
    Parameters
    ----------
    dirpath : str
        Path to a directory
    """
    if dirpath[-1] != '/' and Path(dirpath).is_dir():
        dirpath += '/'
    return dirpath


def parse_cmd_line() -> argparse.Namespace:
    """Parse command line argument.

    Return
    ------
    args : argparse.Namespace
        Object holding validated arguments.
    """
    global DPI, PERFORMANCES_CLASSES
    # Load command line arguments
    args = _get_cmd_line_args()
    # Convert directory paths
    args.benchmark_directory = must_end_with_slash(args.benchmark_directory)
    args.output_path = must_end_with_slash(args.output_path)
    # Set global variables
    DPI = args.dpi
    PERFORMANCES_CLASSES = ALL_PERFORMANCES_CLASSES[args.type]
    # Check that directory exist
    if not os.path.exists(args.benchmark_directory):
        sys.exit(
            f'INPUT ERROR: path to directory to analyse '
            f'"{args.benchmark_directory}" not found!'
            )
    return args


def _get_cmd_line_args() -> argparse.Namespace:
    """Define and parse the command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "benchmark_directory",
        help=(
            "Path to directory "
            "where benchmark was performed by haddock-runner"
            ),
        type=str,
        )
    parser.add_argument(
        "-o",
        "--output_path",
        help="Directory where to write output files.",
        required=False,
        default="Analysis",
        type=str,
        )
    parser.add_argument(
        "-m",
        "--metric",
        help="Performance metric to track.",
        required=False,
        default="irmsd",
        choices=("irmsd", "dockq", ),
        type=str,
        )
    parser.add_argument(
        '-s',
        '--scenario',
        help="Name(s) of a specific scenario(s) to analyze. Can be multiple of them, separated by space. By default, all scenarios will be analysed together.",
        required=False,
        default=None,
        nargs="+",
        type=str,
        )
    parser.add_argument(
        '-t',
        '--type',
        help="Type of analysis to be conducted.",
        required=False,
        default="protein",
        choices=("protein", "peptide", "glycan", ),
        type=str,
        )
    parser.add_argument(
        '-d',
        '--dpi',
        help="DPI of the generated figures.",
        required=False,
        default=400,
        type=int,
        )
    parser.add_argument(
        '--no-capriplots',
        help="Do not generate CAPRI plots",
        action="store_true",
        default=False,
        )
    parser.add_argument(
        '--no-violinplots',
        help="Do not generate violin plots",
        action="store_true",
        default=False,
        )
    parser.add_argument(
        '--no-melquiplots',
        help="Do not generate melqui plots",
        action="store_true",
        default=False,
        )
    parser.add_argument(
        '-n',
        '--no-percentage',
        help="Display number of structures instead of percentages",
        action="store_true",
        default=False,
        )
    parser.add_argument(
        '-q',
        '--quiet',
        help="Silences prints",
        action="store_true",
        default=False,
        )
    args = parser.parse_args()
    return args


def welcome_msg(silent: bool) -> None:
    """Print welcome message

    Parameters
    ----------
    silent : bool
        If true, do not print
    """
    vprint("#" * 80, silent)
    vprint(f"#     Haddock-runner haddock3 analysis script", silent)
    vprint(f"#     Version: {__version__}", silent)
    vprint(f"#     Author:  {__author__}", silent)
    vprint(f"#     Devs:    {', '.join(__dev__)}", silent)
    vprint("#" * 80, silent)


############################
# COMMAND LINE ENTRY POINT #
############################
def maincli() -> None:
    """CLI entry point."""
    args = parse_cmd_line()
    welcome_msg(args.quiet)
    main(
        args.benchmark_directory,
        outputpath=args.output_path,
        scenarios=args.scenario,
        metric=args.metric,
        silent=args.quiet,
        no_percentage=args.no_percentage,
        no_capriplots=args.no_capriplots,
        no_violinplots=args.no_violinplots,
        no_melquiplots=args.no_melquiplots,
        )


if __name__ == "__main__":
    maincli()
