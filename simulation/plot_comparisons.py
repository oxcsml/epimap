import os
import math
from itertools import chain

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_parameter(par_name, results_dir):
    samples = pd.read_csv(os.path.join(results_dir, par_name + ".csv"))
    samples.Date = pd.to_datetime(samples.Date)
    return samples


def load_all(results_dir):

    merged = not os.path.isfile(os.path.join(results_dir, "Rt.csv"))

    df = load_parameter("merged_Rt" if merged else "Rt", results_dir)

    cpred = load_parameter("merged_Cpred" if merged else "Cpred", results_dir)
    cproj = load_parameter("merged_Cproj" if merged else "Cproj", results_dir)
    c = pd.concat([cproj, cpred])

    pexceed = load_parameter("merged_Pexceed" if merged else "Pexceed", results_dir)

    df = pd.merge(df, c, on=["area", "Date", "provenance"])
    df = pd.merge(df, pexceed, on=["area", "Date", "provenance"])

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
    results, par_name, median="_50", upper="_97_5", lower="_2_5", color=None, ax=None
):
    if ax is None:
        fig, ax = plt.subfigures(1, 1)

    if color == None:
        lines = ax.plot(results.Date, results[par_name + median])
        ax.fill_between(
            results.Date,
            results[par_name + lower],
            results[par_name + upper],
            alpha=0.3,
            color=lines[0].get_color(),
        )
    else:
        ax.plot(results.Date, results[par_name + median], color=color)
        ax.fill_between(
            results.Date,
            results[par_name + lower],
            results[par_name + upper],
            alpha=0.3,
            color=color,
        )

    ax.tick_params(axis="x", rotation=45)

    return ax


def plot_actual(results, par_name, color=None, ax=None):
    if ax is None:
        fig, ax = plt.subfigures(1, 1)

    if color == None:
        ax.plot(results.Date, results[par_name])
    else:
        ax.plot(results.Date, results[par_name], color=color)

    ax.tick_params(axis="x", rotation=45)

    return ax


def plot_estimate_area(results, area, *args, **kwargs):
    ax = plot_estimate(results[results.area == area], *args, **kwargs)
    ax.set_title(area)
    return ax


def plot_actual_area(results, area, *args, **kwargs):
    ax = plot_actual(results[results.area == area], *args, **kwargs)
    ax.set_title(area)
    return ax


def plot_comparison(
    actual_results,
    estimates_results,
    par_name,
    title=None,
    axes=None,
    areas=None,
    median="_50",
    upper="_97_5",
    lower="_2_5",
    share_ax=True,
    fix_scale=False,
):
    if areas == None:
        areas = actual_results.area.unique().tolist()

    if axes == None:
        N = len(areas)
        H = math.ceil(math.sqrt(N) / 1.6)
        W = math.ceil(N / H)

        fig, axes = plt.subplots(
            H, W, sharex=share_ax, sharey=share_ax, figsize=(3 * W, 3 * H)
        )

        all_axes = list(chain.from_iterable(axes))
        axes = all_axes[:N]
    else:
        all_axes = axes

    if not isinstance(estimates_results, dict) and estimates_results is not None:
        estimates_results = {"Estimate": estimates_results}

    for area, ax in zip(areas, axes):

        max_y = 0

        if actual_results is not None:
            plot_actual_area(actual_results, area, par_name=par_name, ax=ax)
            max_y = actual_results[actual_results.area == area][par_name].max()

        if estimates_results is not None:
            for estimates_result in estimates_results.values():
                plot_estimate_area(
                    estimates_result,
                    area,
                    par_name=par_name,
                    ax=ax,
                    median=median,
                    upper=upper,
                    lower=lower,
                )
                max_y = max(
                    max_y,
                    estimates_result[estimates_result.area == area][
                        par_name + median
                    ].max(),
                )

        if fix_scale:
            ax.set_ylim([0, 1.3 * max_y])

    for ax in all_axes:
        ax.tick_params(axis="x", rotation=45)

    titles = ["Actual", *list(estimates_results.keys())]
    lines = axes[0].get_lines()

    fig.legend(
        lines,
        titles,
        loc="upper left",
        bbox_to_anchor=(1.0, 0.95),
        bbox_transform=plt.gcf().transFigure,
    )

    if title is not None:
        fig.suptitle(title)

    plt.tight_layout()

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
        fig, axes = plt.subplots(H, W, figsize=(3 * W, 3 * H))

        all_axes = list(chain.from_iterable(axes))
        axes = all_axes[:N]
    else:
        all_axes = axes

    for area, ax in zip(areas, axes):
        plot_actual_area(cases, area, par_name="cases", ax=ax)

    for ax in all_axes:
        ax.tick_params(axis="x", rotation=45)

    plt.tight_layout()

    return fig, axes


if __name__ == "__main__":
    # import os

    # os.chdir("/data/ziz/not-backed-up/mhutchin/Rmap-dev/Rmap")
    # from simulation.plot_comparisons import *

    actuals = load_parameter("Rt", "simulation/latent_epidemic/tehtropolis/sample")
    estimates = {
        "EpiEstim": load_all("simulation_fits/test/epiestim"),
        "SingleArea": load_all("simulation_fits/test/singlearea"),
        "TwoStage": load_all("simulation_fits/test/twostage"),
        "Regional": load_all("simulation_fits/test/regional"),
    }

    cases = load_cases("simulation/latent_epidemic/tehtropolis/sample")

    title = "All"
    fig, axes = plot_comparison(actuals, estimates, title=title, par_name="Rt")
    for ax in axes:
        ax.set_ylim([0, 2])
    fig.savefig(
        f'figures/Rt_{title.lower().replace(" ", "_")}.pdf', bbox_inches="tight"
    )
    fig, axes = plot_comparison(
        cases,
        estimates,
        title=title,
        par_name="C",
        lower="_25",
        upper="_975",
        share_ax=False,
        fix_scale=True,
    )
    fig.savefig(
        f'figures/cases_{title.lower().replace(" ", "_")}.pdf', bbox_inches="tight"
    )
