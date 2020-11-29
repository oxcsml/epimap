from collections import defaultdict
import itertools
import os
from types import SimpleNamespace
import re

import fire
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

import utils


def nicename(run_name):
    name, hps = run_name.split("_", maxsplit=1)
    space = re.search("space_([0-9]*)", hps).group(1)
    time = re.search("time_([0-9]*)", hps).group(1)
    return f"{name[0].upper()}({space},{time})"


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


def align_dates_areas(projections, truth):
    dates = truth.index.intersection(projections.date.unique())
    areas = projections.area.unique()
    return projections.reindex(index=dates), truth[areas].reindex(index=dates)


class Sample:
    def __init__(self, percentiles, truth, pred_key="c_50"):
        self._dates = pd.Index.intersection(percentiles.index, truth.index)
        self.percentiles = percentiles.reindex(index=self._dates)
        self.predictions = self._get_predictions(self.percentiles, pred_key)
        self.truth = truth.reindex(index=self._dates)

    @staticmethod
    def _get_predictions(percentiles, pred_key):
        return (
            projections.reset_index().set_index(["date", "area"])[pred_key].unstack(-1)
        )

    @property
    def delta(self):
        return self.predictions - self.truth

    def logmean_rmse(self, weekly=False):
        return self._logmean_stat(
            err=np.power(self.delta, 2), weekly=weekly, pipefunc=np.sqrt
        )

    def logmean_mae(self, weekly=False):
        return self._logmean_stat(
            err=np.abs(self.delta), weekly=weekly, pipefunc=lambda x: x
        )

    def _logmean_stat(self, err, weekly, pipefunc):
        stat = (
            err.groupby(pd.Grouper(freq="1W"))
            .mean()
            .pipe(pipefunc)
            .reset_index(drop=True)
            .rename(index=lambda x: f"W{x+1}")
            .transpose()
            if weekly
            else err.mean(0)
        )
        return np.log(stat + 1)


def plot_aggregations(all_samples, metrics_dct):
    """
    Args:
        all_samples (dict): dict runs -> area -> weeks -> sample
        metrics_dct (dict): dict str -> func (the metric)
    """
    means = defaultdict(dict)
    for lookback, sample_dct in tqdm(all_samples.items(), desc="Plotting aggregations"):
        fig, axarr = plt.subplots(len(metrics_dct), 1, constrained_layout=True)
        for ax, (mname, metric) in zip(axarr.flatten(), metrics_dct.items()):
            data = pd.concat(utils.map_lowest(metric, sample_dct), axis=1)
            ax.violinplot(data.T, showmeans=True, showextrema=False)
            ax.set_xticks(np.arange(data.columns.shape[0]) + 1)
            ax.set_xticklabels(data.columns, fontsize=7)
            ax.set_ylabel(f"log({mname} + 1)")
            ax.grid(alpha=0.25)
            means[mname][lookback] = data.mean(0)
        ax.set_xlabel("Weeks modelled from 1st June, 2020")
        fig.suptitle(f"History: {lookback} weeks")
        yield fig
        plt.close()
    means = {mname: pd.concat(v, axis=1).T for mname, v in means.items()}

    fig, axarr = plt.subplots(len(means), 1, constrained_layout=True)
    for ax, (mname, df) in zip(axarr.flatten(), means.items()):
        c1 = itertools.cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])
        c2 = itertools.cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])
        for run_name, srs in df.items():
            style = (
                dict(marker="s", color=next(c1), ls="-")
                if run_name.lower().startswith("l")
                else dict(marker="^", color=next(c2), ls=":")
            )
            ax.plot(srs.index, srs, label=run_name, **style)
        ax.set_xticks(weeks)
        ax.set_ylabel(f"log({mname} + 1)")
        ax.grid(alpha=0.25)
    ax.set_xlabel("Weeks modelled from 1st June, 2020")
    axarr.flatten()[0].legend(
        loc="upper left", fontsize=7, ncol=2, title="Lookback History (weeks)"
    )
    fig.suptitle(f"All runs")
    yield fig
    plt.close()


if __name__ == "__main__":
    args = SimpleNamespace()

    def cmd_args(uk_cases, backtests, output):
        """
        Args:
            uk_cases: path to uk_cases.csv
            backtests: path to backtests folder
            output: folder to save outputs
        """
        args.uk_cases_pth = uk_cases
        args.backtests_dir = backtests
        args.outputs_dir = output

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

    uk_cases = load_uk_cases(args.uk_cases_pth)
    all_samples = defaultdict(dict)
    for folder, week in itertools.product(run_names, weeks):
        projections = read_csv(
            os.path.join(
                args.backtests_dir,
                folder,
                f"map_start_2020-06-01_weeks_{week}",
                "merged_Cproj.csv",
            )
        )

        all_samples[week][nicename(folder)] = Sample(
            *align_dates_areas(projections, uk_cases), pred_key="c_50"
        )

    metrics_dct = {
        "RMSE": lambda x: x.logmean_rmse(weekly=False),
        "MAE": lambda x: x.logmean_mae(weekly=False),
    }

    deck = utils.PdfDeck()
    for fig in plot_aggregations(all_samples, metrics_dct):
        deck.add_figure(fig)
    deck.figs.insert(0, deck.figs.pop(-1))  # put the average over areas one at the front
    print("Writing pdf...")
    deck.make(os.path.join(args.outputs_dir, "backtest_evaluation_aggregations.pdf"))

    print("Making table...")
    weekly_rmse = utils.map_lowest(
        lambda x: x.logmean_rmse(weekly=True).mean(0), all_samples
    )
    weekly_mae = utils.map_lowest(
        lambda x: x.logmean_mae(weekly=True).mean(0), all_samples
    )

    out_table = pd.concat(
        {
            "RMSE": utils.stack_table(weekly_rmse).unstack(-1),
            "MAE": utils.stack_table(weekly_mae).unstack(-1),
        },
        axis=1,
    ).rename_axis(index=["lookback", "counts(time,space)"])
    out_table.to_csv(os.path.join(args.outputs_dir, "table.csv"))
