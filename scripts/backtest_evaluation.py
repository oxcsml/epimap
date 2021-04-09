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

# need a better solution for this in the long run..
import utils

# import scripts.utils as utils


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
    def __init__(self, predictions, truth):
        self._dates = pd.Index.intersection(predictions.index, truth.index)
        self.predictions = predictions.reindex(index=self._dates)
        self.truth = truth.reindex(index=self._dates)

        self.delta = self.predictions - self.truth
        # Alternative option?
        # self.log_delta = np.log(self.predictions + 1) - np.log(self.truth + 1)

    def logmean_rmse(self):
        return self._logmean_stat(err=np.power(self.delta, 2), pipefunc=np.sqrt)

    def logmean_mae(self):
        return self._logmean_stat(err=np.abs(self.delta), pipefunc=None)

    def _logmean_stat(self, err, pipefunc=None):
        grouped = err.groupby(pd.Grouper(freq="1W")).mean()
        grouped = pipefunc(grouped) if pipefunc is not None else grouped
        weekly = (
            grouped.reset_index(drop=True).rename(index=lambda x: f"W{x+1}").transpose()
        )
        means = pipefunc(err.mean(0)) if pipefunc is not None else err.mean(0)
        means.name = "Mean"
        return np.log(pd.concat([weekly, means], axis=1) + 1)


def bar_chart(ax, df, colors):
    c1 = itertools.cycle(colors)
    n_run_names = df.index.get_level_values("model").nunique()
    dates = df.index.get_level_values("start_date").unique().values
    width = 1 / (n_run_names + 1)

    for run_name_idx, (run_name, srs) in enumerate(df.groupby(level="model")):
        srs = srs.droplevel("model").squeeze()
        style = dict(color=next(c1), width=width)
        offset = width * run_name_idx - 0.5 + width
        ax.bar(np.arange(len(srs)) + offset, srs, label=run_name, **style)
    ax.set_xticks(np.arange(len(dates)))
    ax.set_xticklabels(dates)
    ax.grid(alpha=0.25)


def violinplot(ax, data):
    ax.violinplot(data, showmeans=True, showextrema=False)
    ax.set_xticks(np.arange(data.columns.shape[0]) + 1)
    ax.set_xticklabels(data.columns, fontsize=7)
    ax.grid(alpha=0.25)


def scatter(ax, x, y):
    ax.scatter(x, y, marker=".", alpha=0.25)
    ax.grid(alpha=0.25)


def make_plots(df, bar_colors):
    nmetrics = df.columns.get_level_values("metric").nunique()

    # the bar charts
    area_means = df.mean(0).to_frame()
    for freq, stats in area_means.groupby(level="freq"):
        fig, axarr = plt.subplots(nmetrics, 1, constrained_layout=True)
        for ax, (metric, results) in zip(axarr.flat, stats.groupby(level="metric")):
            bar_chart(ax, results.droplevel(["freq", "metric"]), colors=bar_colors)
            ax.set_ylabel(f"{metric} of log(# + 1)")
        ax.set_xlabel("Start date of 15 week modelled period.")
        axarr.flatten()[0].legend(
            loc="upper left", fontsize=7, ncol=2, title="Models", framealpha=0.5
        )
        fig.suptitle(f"All runs: {freq}")
        yield fig
        plt.close()

    # scatter plots should go here
    comparator = "stage1"
    nomean = df.transpose().drop("Mean", level="freq")
    levels = ["metric", "start_date", "freq"]
    for (metric, start_date, freq), res in nomean.groupby(level=levels):
        stats = res.droplevel(levels)
        xaxis = stats.loc[comparator]
        stats = stats.drop(comparator)
        fig, axarr = plt.subplots(3, 3, constrained_layout=True, sharex=True)
        for ax, (method, results) in zip(axarr.flat, stats.iterrows()):
            scatter(ax, xaxis, results)
            ax.set_xlabel(comparator)
            ax.set_ylabel(method)
        fig.suptitle(f"{metric}, start date: {start_date}, {freq}")
        yield fig
        plt.close()

    # violin plots
    mow = df.transpose().query("freq=='Mean'").droplevel("freq")
    for start_date, stats in mow.groupby(level="start_date"):
        fig, axarr = plt.subplots(nmetrics, 1, constrained_layout=True)
        for ax, (metric, results) in zip(axarr.flat, stats.groupby(level="metric")):
            data = results.droplevel(["start_date", "metric"]).transpose()
            violinplot(ax, data)
            ax.set_ylabel(f"{metric} of log(# + 1)")
        ax.set_xlabel("Model")
        forecast_start = (
            pd.Timestamp(start_date) + pd.Timedelta(7, unit="W")
        ).strftime("%Y-%m-%d")
        fig.suptitle(
            (
                "Modelled period: "
                f"{start_date} + 15 weeks.\nPredictions period: {forecast_start} + 3 weeks."
                "\nViolins across regional errors."
            )
        )
        yield fig
        plt.close()


if __name__ == "__main__":

    class Args:
        def populate(self, uk_cases, backtests, output):
            """
            Args:
                uk_cases: path to uk_cases.csv
                backtests: path to backtests folder
                output: folder to save outputs
            """
            self.uk_cases_pth = uk_cases
            self.backtests_dir = backtests
            self.outputs_dir = output

    args = Args()
    fire.Fire(args.populate)

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
        for name in run_names
        if name.startswith("space")
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

    print("Organising data...")
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
            aligned_projections.iloc[0] = uk_cases[aligned_projections.columns].loc[
                aligned_projections.index[0]
            ]
            aligned_projections.ffill(inplace=True)

        all_samples[start_date][folder] = Sample(
            aligned_projections,
            aligned_true_cases,
        )

    data = {
        "RMSE": utils.map_lowest(lambda x: x.logmean_rmse(), all_samples),
        "MAE": utils.map_lowest(lambda x: x.logmean_mae(), all_samples),
    }

    data = utils.collapse(data, axis=1, names=["metric", "start_date", "model", "freq"])

    print("Making charts...")
    deck = utils.PdfDeck()
    for fig in make_plots(data, bar_colors=colors):
        deck.add_figure(fig)

    deck.make(
        os.path.join(
            args.outputs_dir,
            f"backtest_evaluation_aggregations.pdf",
        )
    )

    print("Making table...")
    out_table = (
        data.mean(0)
        .unstack(0)
        .unstack(-1)
        .rename_axis(["", ""], axis=1)
        .rename_axis(["lookback", "counts(time,space)"], axis=0)
    )

    out_table.to_csv(os.path.join(args.outputs_dir, "table.csv"))
