# %%
import os
import math
from itertools import chain
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def lighten_color(color, amount=1.0):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys

    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


def load_parameter(par_name, results_dir):
    samples = pd.read_csv(os.path.join(results_dir, par_name + ".csv"))
    samples.Date = pd.to_datetime(samples.Date)
    return samples


def load_all(results_dir):

    merged = not os.path.isfile(os.path.join(results_dir, "Rt.csv"))

    df = load_parameter("merged_Rt" if merged else "Rt", results_dir)

    # cpred = load_parameter("merged_Cpred" if merged else "Cpred", results_dir)
    # cproj = load_parameter("merged_Cproj" if merged else "Cproj", results_dir)
    # c = pd.concat([cproj, cpred])
    # df = pd.merge(df, c, on=["area", "Date", "provenance"])

    bpred = load_parameter("merged_Bpred" if merged else "Bpred", results_dir)
    bproj = load_parameter("merged_Bproj" if merged else "Bproj", results_dir)
    b = pd.concat([bproj, bpred])
    df = pd.merge(df, b, on=["area", "Date", "provenance"])

    # xpred = load_parameter("merged_Xpred" if merged else "Xpred", results_dir)
    # xproj = load_parameter("merged_Xproj" if merged else "Xproj", results_dir)
    # x = pd.concat([xproj, xpred])
    # df = pd.merge(df, x, on=["area", "Date", "provenance"])

    # pexceed = load_parameter("merged_Pexceed" if merged else "Pexceed", results_dir)
    # df = pd.merge(df, pexceed, on=["area", "Date", "provenance"])

    data_cols = df.columns.drop(["area", "Date", "provenance"])
    df[data_cols] = df[data_cols].replace({"  NA": "NaN"}).astype(np.float32)

    return df


def load_cases(data_dir):
    cases = pd.read_csv(os.path.join(data_dir, "cases.csv"))
    cases = cases.melt(id_vars="area", var_name="Date", value_name="cases").sort_values(
        by=["area", "Date"]
    )
    cases.Date = pd.to_datetime(cases.Date)
    cases = cases.rename(columns={"cases": "C"})
    return cases


def plot_estimate(
    results,
    par_name,
    provenance=None,
    median="_50",
    upper1="_97_5",
    lower1="_2_5",
    upper2="_97_5",
    lower2="_2_5",
    color=None,
    ax=None,
):
    if ax is None:
        fig, ax = plt.subfigures(1, 1)

    if provenance is not None:
        results = results[results.provenance == provenance]

    if provenance == "projected":
        linestyle = "dotted"
    else:
        linestyle = "solid"

    if color is None:
        line = ax.plot(
            results.Date,
            results[par_name + median],
            linewidth=1.0,
            linestyle=linestyle,
            zorder=1,
        )
        color = (line[0].get_color(),)
        fill1 = ax.fill_between(
            results.Date,
            results[par_name + lower1],
            results[par_name + upper1],
            alpha=0.2 if provenance == "projected" else 0.4,
            color=color,
            zorder=2,
        )
        fill2 = ax.fill_between(
            results.Date,
            results[par_name + lower2],
            results[par_name + upper2],
            alpha=0.2 if provenance == "projected" else 0.4,
            color=color,
            zorder=2,
        )
    else:
        line = ax.plot(
            results.Date,
            results[par_name + median],
            color=color,
            linewidth=1.0,
            linestyle=linestyle,
            zorder=1,
        )
        fill1 = ax.fill_between(
            results.Date,
            results[par_name + lower1],
            results[par_name + upper1],
            alpha=0.4 if provenance == "projected" else 0.4,
            color=color,
            zorder=2,
            # hatch="///" if provenance == "projected" else None,
            # facecolor=lighten_color(color, 0.5),
            # edgecolor=lighten_color(color, 1.6),
            # linewidth=0.0,
        )
        fill2 = ax.fill_between(
            results.Date,
            results[par_name + lower2],
            results[par_name + upper2],
            alpha=0.2 if provenance == "projected" else 0.4,
            color=color,
            zorder=2,
            # hatch="///" if provenance == "projected" else None,
            # facecolor=lighten_color(color, 0.5),
            # edgecolor=lighten_color(color, 1.6),
            # linewidth=0.0,
        )

    ax.tick_params(axis="x", rotation=45)

    return (line, fill1, fill2), ax


def plot_actual(results, par_name, color=None, scatter=False, ax=None):
    if ax is None:
        fig, ax = plt.subfigures(1, 1)

    # if color == None:
    #     # ax.plot(results.Date, results[par_name], zorder=1)
    #     ax.scatter(results.Date, results[par_name], zorder=1)
    # else:
    if scatter:
        ax.scatter(
            results.Date,
            results[par_name],
            marker="+",
            s=4,
            color=color,
            zorder=1,
        )
    else:
        ax.plot(results.Date, results[par_name], color=color, zorder=1)

    ax.tick_params(axis="x", rotation=45)

    return ax


def plot_estimate_area(results, area, *args, **kwargs):
    handles, ax = plot_estimate(results[results.area == area], *args, **kwargs)
    # ax.set_title(area)
    return handles, ax


def plot_actual_area(results, area, *args, **kwargs):
    ax = plot_actual(results[results.area == area], *args, **kwargs)
    # ax.set_title(area)
    return ax


def plot_comparison(
    actual_results,
    estimates_results,
    par_name,
    title=None,
    axes=None,
    areas=None,
    median="_50",
    upper1="_75",
    lower1="_25",
    upper2="_97_5",
    lower2="_2_5",
    scatter_actual=False,
    sharex=True,
    sharey=True,
    fix_scale=False,
    estimate_colors=["tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"],
):
    if areas == None:
        areas = actual_results.area.unique().tolist()

    if axes == None:
        N = len(areas)
        W = len(areas)
        H = len(estimates_results)

        fig, axes = plt.subplots(
            H,
            W,
            sharex=sharex,
            sharey=sharey,
            figsize=(2.5 * W, 1.5 * H),
            squeeze=False,
        )

        all_axes = list(chain.from_iterable(axes))
    else:
        all_axes = list(chain.from_iterable(axes))

    if not isinstance(estimates_results, dict) and estimates_results is not None:
        estimates_results = {"Estimate": estimates_results}

    # max_y = {i: 0.0 for i in range(W)}
    max_y = defaultdict(float)

    handles = []
    labels = []

    for i, (row, (name, estimate)) in enumerate(zip(axes, estimates_results.items())):
        for j, (area, ax) in enumerate(zip(areas, row)):

            if actual_results is not None:
                plot_actual_area(
                    actual_results,
                    area,
                    scatter=scatter_actual,
                    par_name=par_name,
                    ax=ax,
                )
                max_y[j] = max(
                    max_y[j],
                    actual_results[actual_results.area == area][par_name].max(),
                )

            infr_handles, ax = plot_estimate_area(
                estimate,
                area,
                par_name=par_name,
                ax=ax,
                median=median,
                upper1=upper1,
                lower1=lower1,
                upper2=upper2,
                lower2=lower2,
                color=estimate_colors[i],
                provenance="inferred",
            )
            proj_handles, ax = plot_estimate_area(
                estimate,
                area,
                par_name=par_name,
                ax=ax,
                median=median,
                upper1=upper1,
                lower1=lower1,
                upper2=upper2,
                lower2=lower2,
                color=estimate_colors[i],
                provenance="projected",
            )
            ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
            max_y[j] = max(
                max_y[j],
                estimate[estimate.area == area][par_name + median].max(),
            )

            if i == 0:
                ax.set_title(f"{area}")

            if j == 0:
                handles.append(infr_handles)
                handles.append(proj_handles)
                labels += [f"{name} infered", f"{name} projected"]

    for i, row in enumerate(axes):
        for j, ax in enumerate(row):
            ax.tick_params(axis="x", rotation=45)
            if fix_scale:
                ax.set_ylim([0, 1.3 * max_y[j]])

    titles = ["Actual", *list(estimates_results.keys())]
    # lines = axes[0].get_lines()

    # fig.legend(
    #     lines,
    #     titles,
    #     loc="upper left",
    #     bbox_to_anchor=(1.0, 0.95),
    #     bbox_transform=plt.gcf().transFigure,
    # )

    if title is not None:
        fig.suptitle(title)

    # print(handles)
    # print([h[1] for h in handles])
    # print(labels)
    plt.tight_layout()

    # for h in handles:
    #     print(h)

    plt.legend(
        [(h[0][0], h[1], h[2]) for h in handles],
        labels,
        loc="upper left",
        bbox_to_anchor=(1.0, 0.95),
        bbox_transform=plt.gcf().transFigure,
    )

    return fig, axes


def plot_cases(
    cases,
    axes=None,
    areas=None,
):
    if areas == None:
        areas = cases.area.unique().tolist()

    if axes == None:
        N = len(areas)
        H = math.ceil(math.sqrt(N) / 1.6)
        W = math.ceil(N / H)

        # fig, axes = plt.subplots(H, W, sharex=True, sharey=True, figsize=(3 * W, 3 * H))
        fig, axes = plt.subplots(H, W, figsize=(3 * W, 3 * H), squeeze=False)

        all_axes = list(chain.from_iterable(axes))
        # axes = all_axes[:N]
    else:
        all_axes = list(chain.from_iterable(axes))

    for row in axes:
        for area, ax in zip(areas, row):
            plot_actual_area(cases, area, par_name="cases", ax=ax)

    for ax in all_axes:
        ax.tick_params(axis="x", rotation=45)

    plt.tight_layout()

    return fig, axes


if __name__ == "__main__":
    import os

    os.chdir("/data/ziz/not-backed-up/mhutchin/Rmap-dev/Rmap")
    # from simulation.plot_comparisons import *
    sample_folder = "simulation/latent_epidemic/tehtropolis/sample_rss_paper_1/"
    results_folder = "simulation_fits/sample_rss_paper_1_longer/"
    actuals = load_parameter("Rt", sample_folder)
    estimates = {
        "EpiMap (single area)": load_all(results_folder + "singlearea"),
        # "TwoStage": load_all(results_folder + "twostage"),
        "EpiMap (spatial: 0.05, time: 200)": load_all(results_folder + "regional"),
        "EpiMap (spatial: 0.5, time: 200)": load_all(
            "simulation_fits/sample_rss_paper_1_longer_space/" + "regional"
        ),
        "EpiNow2": load_all(results_folder + "epinow2"),
        "EpiEstim": load_all(results_folder + "epiestim"),
    }

    cases = load_cases(sample_folder)

    title = "All"

    mapping = {
        "Tehtropolis": "Oxford",
        "Brynshire": "Cherwell",
        "Bobbingdon": "West Oxfordshire",
        "Hutchintown": "South Oxfordshire",
        "Shehland": "Buckinghamshire",
    }

    cases.replace({"area": mapping}, inplace=True)
    actuals.replace({"area": mapping}, inplace=True)
    for estimate in estimates.values():
        estimate.replace({"area": mapping}, inplace=True)

    sweep_colors = list(matplotlib.cm.get_cmap("plasma", 9).colors)[-8:-1]

    estimate_colors = [
        "black",
        sweep_colors[2],
        sweep_colors[5],
        "forestgreen",
        "crimson",
    ]
    # %%
    fig, axes = plot_comparison(
        actuals,
        estimates,
        par_name="Rt",
        estimate_colors=estimate_colors,
    )
    for ax in chain.from_iterable(axes):
        ax.set_ylim([0, 3])
        ax.yaxis.get_major_locator().set_params(integer=True)
        # ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    # plt.locator_params(nbins=10)

    fig.savefig(
        f'figures/Rt_{title.lower().replace(" ", "_")}.pdf', bbox_inches="tight"
    )
    fig, axes = plot_comparison(
        cases,
        estimates,
        par_name="C",
        lower1="_25",
        upper1="_75",
        lower2="_025",
        upper2="_975",
        sharey=False,
        fix_scale=True,
        scatter_actual=True,
        estimate_colors=estimate_colors,
    )
    # for ax in chain.from_iterable(axes):
    # ax.set_yscale("log")
    fig.savefig(
        f'figures/cases_{title.lower().replace(" ", "_")}.pdf', bbox_inches="tight"
    )
    # fig, axes = plot_comparison(
    #     cases,
    #     estimates,
    #     par_name="B",
    #     lower="_25",
    #     upper="_975",
    #     sharey=False,
    #     fix_scale=True,
    # )
    # fig.savefig(
    #     f'figures/cases_no_weekly_{title.lower().replace(" ", "_")}.pdf',
    #     bbox_inches="tight",
    # )
