from collections import defaultdict
import itertools
from functools import partial
import os
import pickle
from types import SimpleNamespace

import fire
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

import utils

"""
Note: 
    Rather than making one sample per area, it would be better to just 
    make one sample per run and take areas out for plotting...


    If we use this code again we can make that change.
    It would just mean replacing the call to group samples with
    Sample(predictions=predictions, truth=uk_cases, ...)

    The script will run a _lot_ faster like this!
"""


def read_csv(pth):
    df = pd.read_csv(pth).rename(columns=str.lower)
    df.date = pd.to_datetime(df.date, format="%Y-%m-%d")
    return df


def load_uk_cases(pth):
    uk_cases = pd.read_csv(pth).set_index("Area name").drop("Country", 1).transpose()
    uk_cases.index = pd.to_datetime(uk_cases.index)
    uk_cases.index.name = "date"
    uk_cases.columns.name = None
    return uk_cases


def stack_table(dct, axis=0, stack_func=None):
    "recursively concat dict of dict of ... of dataframes"
    stacker = stack_func or partial(pd.concat, axis=axis)
    return stacker(
        {
            k: stack_table(v)
            if isinstance(next(iter(v.values())), dict)
            else stacker(v)
            for k, v in dct.items()
        }
    )


class Sample:
    def __init__(self, percentiles, truth, pred_key="c_50"):
        self._dates = pd.Index.intersection(percentiles.index, truth.index)
        self.percentiles = percentiles.reindex(index=self._dates)
        self.predictions = self.percentiles[pred_key]
        self.truth = truth.reindex(index=self._dates)


def group_samples(counts, preds, index, by, return_missing=False):
    grouped_preds = preds.set_index(index).groupby(by)
    all_samples = dict()
    missing_from_predictions = list()
    for area in counts.columns:
        try:
            pred = grouped_preds.get_group(area)
        except KeyError:
            missing_from_predictions.append(area)
            continue
        else:
            all_samples[area] = Sample(pred.drop([by], axis=1), counts[area].squeeze(),)
    return (all_samples, missing_from_predictions) if return_missing else all_samples


def rmse(predictions, truth):
    return np.sqrt(np.mean(np.power(predictions - truth, 2)))


def mae(predictions, truth):
    return np.mean(np.abs(predictions - truth))


def get_weekly_rmse(sample):
    delta = np.power(sample.predictions - sample.truth, 2)
    return np.sqrt(delta.groupby(pd.Grouper(freq="1W")).mean()).reset_index(drop=True)


def get_weekly_mae(sample):
    delta = np.abs(sample.predictions - sample.truth)
    return np.sqrt(delta.groupby(pd.Grouper(freq="1W")).mean()).reset_index(drop=True)


def apply_on_each(metric, samples):
    return [metric(s.predictions, s.truth) for s in samples.values()]


def plot_all_areas(all_samples, metrics_dct):
    """
    Args:
        all_samples (dict): dict area -> runs -> weeks -> sample
        metrics_dct (dict): dict str -> func (the metric)
    """
    for area, sample_dct in tqdm(all_samples.items(), desc="Plotting all areas"):
        fig, axarr = plt.subplots(
            len(metrics_dct),
            1,
            figsize=(10, 3 * len(metrics_dct)),
            constrained_layout=True,
            sharex=True,
        )
        for ax, mname in zip(axarr, metrics_dct):
            for run, samples in sample_dct.items():
                ax.plot(
                    list(samples.keys()),
                    apply_on_each(metrics_dct[mname], samples),
                    label=run,
                )
            ax.legend(fontsize=7, ncol=2)
            ax.grid(alpha=0.25)
            ax.set_ylabel(mname)

        axarr[-1].set_xlabel("Weeks modelled from 1st June, 2020")
        fig.suptitle(area)
        yield fig
        plt.close()


def plot_aggregations(all_samples, metrics_dct):
    """
    Args:
        all_samples (dict): dict runs -> area -> weeks -> sample
        metrics_dct (dict): dict str -> func (the metric)
    """
    logmeans = defaultdict(dict)
    for run_name, sample_dct in tqdm(all_samples.items(), desc="Plotting aggregations"):
        weeks = list(next(iter(sample_dct.values())).keys())
        fig, axarr = plt.subplots(len(metrics_dct), 1, constrained_layout=True)
        for ax, (mname, metric) in zip(axarr.flatten(), metrics_dct.items()):
            # data = np.row_stack([apply_on_each(metric, s) for s in sample_dct.values()])
            # transdata = np.log(data + 1)
            # logmeans[metric][run_name] = transdata.mean(0)
            ax.violinplot(transdata, showmeans=True, showextrema=False, positions=weeks)
            ax.set_xticks(weeks)
            ax.set_ylabel(f"log({mname} + 1)")
            ax.grid(alpha=0.25)
        ax.set_xlabel("Weeks modelled from 1st June, 2020")
        fig.suptitle(f"{run_name}")
        yield fig
        plt.close()

    fig, axarr = plt.subplots(len(metrics_dct), 1, constrained_layout=True)
    for ax, (mname, metric) in zip(axarr.flatten(), metrics_dct.items()):
        c1 = itertools.cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])
        c2 = itertools.cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])
        for run_name, data in logmeans[metric].items():
            style = (
                dict(marker="s", color=next(c1), ls="-")
                if run_name.startswith("latent")
                else dict(marker="^", color=next(c2), ls=":")
            )
            ax.plot(weeks, data, label=run_name, **style)
            ax.set_xticks(weeks)
        ax.set_ylabel(f"log({mname} + 1)")
        ax.grid(alpha=0.25)
    ax.set_xlabel("Weeks modelled from 1st June, 2020")
    axarr.flatten()[0].legend(loc="upper left", fontsize=7, ncol=2)
    fig.suptitle(f"All runs")
    yield fig
    plt.close()


if __name__ == "__main__":
    args = SimpleNamespace()

    def cmd_args(
        uk_cases,
        backtests,
        output,
        plot_all_areas=False,
        plot_aggregations=True,
        table=True,
        pkl_load=None,
        pkl_dump=None,
    ):
        """
        Args:
            uk_cases: path to uk_cases.csv
            backtests: path to backtests folder
            output: folder to save outputs
        """
        args.uk_cases_pth = uk_cases
        args.backtests_dir = backtests
        args.outputs_dir = output
        args.table = table
        args.pkl_load = pkl_load
        args.pkl_dump = pkl_dump
        args.plot_all_areas = plot_all_areas
        args.plot_aggregations = plot_aggregations

    fire.Fire(cmd_args)

    # todo: just infer this from the filenames
    weeks = [4, 8, 12, 16]
    run_names = [
        "latents_space_1_time_60",
        "latents_space_3_time_30",
        "latents_space_3_time_60",
        "latents_space_3_time_120",
        "latents_space_5_time_60",
        "reports_space_3_time_30",
        "reports_space_1_time_60",
        "reports_space_3_time_60",
        "reports_space_3_time_120",
        "reports_space_5_time_60",
    ]

    # could thread...?
    if args.pkl_load is not None:
        with open(args.pkl_load, "rb") as f:
            all_samples = pickle.load(f)
    else:
        uk_cases = load_uk_cases(args.uk_cases_pth)
        all_samples = defaultdict(lambda: defaultdict(dict))
        for folder, week in itertools.product(run_names, weeks):
            projections = read_csv(
                os.path.join(
                    args.backtests_dir,
                    folder,
                    f"map_start_2020-06-01_weeks_{week}",
                    "merged_Cproj.csv",
                )
            )

            all_samples[folder][week] = Sample(
                    percentiles=projections,
                    truth=uk_cases,
                    pred_key="c_50"
                )

    if args.pkl_dump is not None:
        with open(args.pkl_dump, "wb") as f:
            pickle.dump(all_samples, f, pickle.HIGHEST_PROTOCOL)

    metrics_dct = {"RMSE": rmse, "MAE": mae}

    if args.table:
        print("Making table...")
        weekly_rmse = stack_table(utils.map_lowest(get_weekly_rmse, all_samples))
        weekly_mae = stack_table(utils.map_lowest(get_weekly_mae, all_samples))

        def logmean_pretty(long_table):
            return (
                np.log(long_table + 1)
                .unstack(level=1)
                .mean(axis=1)
                .unstack(level=-1)
                .rename(columns=lambda x: f"W{x+1}")
            )
        out_table = pd.concat(
            {"RMSE": logmean_pretty(weekly_rmse), "MAE": logmean_pretty(weekly_mae)}, axis=1
        )
        out_table.to_csv(os.path.join(args.outputs_dir, "table.csv"))

