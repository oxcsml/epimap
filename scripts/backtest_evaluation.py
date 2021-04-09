from collections import defaultdict
import itertools
import os
from types import SimpleNamespace
import re

import fire
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

import utils

# import scripts.utils as utils


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


def _get_predictions(percentiles, pred_key):
    return percentiles.reset_index().set_index(["date", "area"])[pred_key].unstack(-1)


def align_dates_areas(projections, truth, pred_key="c_50"):
    dates = truth.index.intersection(pd.DatetimeIndex(projections.date.unique()))
    areas = projections.area.unique()
    predictions = _get_predictions(projections, pred_key)
    return predictions[areas].reindex(index=dates), truth[areas].reindex(index=dates)


class Sample:
    def __init__(self, percentiles, truth):  # pred_key="c_50"
        self._dates = pd.Index.intersection(percentiles.index, truth.index)
        self.predictions = percentiles.reindex(index=self._dates)
        # self.predictions = self._get_predictions(self.percentiles, pred_key)
        self.truth = truth.reindex(index=self._dates)

    # @staticmethod
    # def _get_predictions(percentiles, pred_key):
    #     return (
    #         projections.reset_index().set_index(["date", "area"])[pred_key].unstack(-1) # SZ: doesn't depend on percentiles?
    #     )

    @property
    def delta(self):
        return (
            self.predictions - self.truth
        )  # np.log(self.predictions + 1) - np.log(self.truth + 1) #

    def logmean_rmse(self, weekly=False):
        return self._logmean_stat(
            err=np.power(self.delta, 2), weekly=weekly, pipefunc=np.sqrt
        )

    def logmean_mae(self, weekly=False):
        return self._logmean_stat(
            err=np.abs(self.delta), weekly=weekly, pipefunc=lambda x: x
        )

    def _logmean_stat(
        self, err, weekly, pipefunc
    ):  # SZ: no square root if weekly = False?
        stat = (
            err.groupby(pd.Grouper(freq="1W"))
            .mean()
            .pipe(pipefunc)
            .reset_index(drop=True)
            .rename(index=lambda x: f"W{x+1}")
            .transpose()
            if weekly
            else err.mean(0).pipe(pipefunc)
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
            ax.violinplot(data, showmeans=True, showextrema=False)
            ax.set_xticks(np.arange(data.columns.shape[0]) + 1)
            ax.set_xticklabels(data.columns, fontsize=7)
            # ax.set_ylabel(f"log({mname} + 1)")
            ax.set_ylabel(f"{mname} of log(# + 1)")
            ax.grid(alpha=0.25)
            means[mname][lookback] = data.mean(0)
        ax.set_xlabel("Model")
        predictions_start = str(sample_dct["stage1"]._dates[0].date())
        fig.suptitle(
            f"Modelled period: {lookback} + 15 weeks. \nPredictions period: {predictions_start} + 3 weeks. \nViolins across regional errors."
        )
        yield fig
        plt.close()
    means = {mname: pd.concat(v, axis=1).T for mname, v in means.items()}

    fig, axarr = plt.subplots(len(means), 1, constrained_layout=True)
    for ax, (mname, df) in zip(axarr.flatten(), means.items()):
        # c1 = itertools.cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])
        c1 = itertools.cycle(colors)
        # c2 = itertools.cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])

        for run_name_idx, (run_name, srs) in enumerate(df.items()):
            if line_plot:
                style = (
                    dict(
                        marker="o", color=next(c1), ls="-", linewidth=0.5, markersize=2
                    )
                    # if run_name.lower().startswith("l")
                    # else dict(marker="^", color=next(c2), ls=":")
                )
                ax.plot(srs.index, srs, label=run_name, **style)
            else:
                style = (
                    dict(color=next(c1), width=width)
                    # if run_name.lower().startswith("l")
                    # else dict(marker="^", color=next(c2), ls=":")
                )
                offset = width * run_name_idx - 0.5 + width
                ax.bar(np.arange(len(srs.index)) + offset, srs, label=run_name, **style)
        ax.set_xticks(np.arange(len(srs.index)))
        ax.set_xticklabels(srs.index)
        # ax.set_ylabel(f"log({mname} + 1)")
        ax.set_ylabel(f"{mname} of log(# + 1)")
        ax.grid(alpha=0.25)
    ax.set_xlabel("Start date of 15 week modelled period.")
    axarr.flatten()[0].legend(
        loc="upper left", fontsize=7, ncol=2, title="Models", framealpha=0.5
    )
    fig.suptitle(f"All runs")
    yield fig
    plt.close()


if __name__ == "__main__":

    args = SimpleNamespace()

    def cmd_args(uk_cases, backtests, output, weekly=False):
        """
        Args:
            uk_cases: path to uk_cases.csv
            backtests: path to backtests folder
            output: folder to save outputs
        """
        args.uk_cases_pth = uk_cases
        args.backtests_dir = backtests
        args.outputs_dir = output
        args.weekly = weekly

    fire.Fire(cmd_args)

    # todo: just infer this from the filenames
    start_dates = [
        "2020-08-24",
        "2020-09-07",
        "2020-09-21",
        "2020-10-05",
        "2020-10-19",
    ]

    run_names = [
        "stage1",
        "space_0",
        "space_0.01",
        "space_0.05",
        "space_0.1",
        "space_0.2",
        "space_0.5",
        "space_1.0",
        "zeros",
        "last_case_count",
    ]
    stub = "start_{start_date}_weeks_{weeks}"
    path_dct = {
        name: os.path.join(
            args.backtests_dir,
            name,
            stub,
            "regional",
            "merged_Cproj.csv",
        )
        for name in run_names if name.startswith("space")
    }
    path_dct.update(
        {
            "stage1": os.path.join(
                args.backtests_dir,
                "space_0",
                stub,
                "singlearea",
                "Cproj.csv",
            ),
            "zeros": os.path.join(
                args.backtests_dir,
                "space_0",
                stub,
                "singlearea",
                "Cproj.csv",
            ),
            # any df in right format is ok here
            "last_case_count": os.path.join(
                args.backtests_dir,
                "space_0",
                stub,
                "singlearea",
                "Cproj.csv",
            ),
        }
    )
    colors = (
        ["black"]
        + list(matplotlib.cm.get_cmap("viridis", 7).colors)
        + ["crimson", "dodgerblue"]
    )
    line_plot = False
    width = 1 / (len(run_names) + 1)
    pred_key = "c_50"

    uk_cases = load_uk_cases(args.uk_cases_pth)
    all_samples = defaultdict(dict)
    for folder, start_date in itertools.product(run_names, start_dates):
        projections = read_csv(path_dct[folder].format(start_date=start_date, weeks=15))

        if folder == "zeros":
            projections[pred_key] = 0.0

        aligned_projections, aligned_true_cases = align_dates_areas(
            projections, uk_cases, pred_key=pred_key
        )

        if folder == "last_case_count":
            aligned_projections *= np.nan
            aligned_projections.iloc[0] = uk_cases[aligned_projections.columns].loc[aligned_projections.index[0]]
            aligned_projections.ffill(inplace=True)

        all_samples[start_date][folder] = Sample(
            aligned_projections,
            aligned_true_cases,
        )

    # Weekly: split one week per chart so you can see it
    #  maybe also just do this automatically, then add the extra charts in the deck
    #  and take out the weekly flag

    metrics_dct = {
        "RMSE": lambda x: x.logmean_rmse(weekly=args.weekly),
        "MAE": lambda x: x.logmean_mae(weekly=args.weekly),
    }

    deck = utils.PdfDeck()
    for fig in plot_aggregations(all_samples, metrics_dct):
        deck.add_figure(fig)
    deck.figs.insert(
        0, deck.figs.pop(-1)
    )  # put the average over areas one at the front
    print("Writing pdf...")
    deck.make(
        os.path.join(
            args.outputs_dir,
            f"backtest_evaluation_aggregations{'_weekly' if args.weekly else ''}.pdf",
        )
    )

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
