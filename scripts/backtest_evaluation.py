from collections import defaultdict
import itertools
import os
import pickle
from types import SimpleNamespace

import fire
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

import utils


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


class Sample:
    def __init__(self, percentiles, truth, pred_key="c_50", slicer=None):
        self._dates = pd.Index.intersection(percentiles.index, truth.index)
        if slicer is not None:
            self._dates = self._dates[slicer]
        self.percentiles = percentiles.reindex(index=self._dates)
        self.predictions = self.percentiles[pred_key]
        self.truth = truth.reindex(index=self._dates)


def group_samples(counts, preds, index, by, return_missing=False, slicer=None):
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
            all_samples[area] = Sample(
                pred.drop([by], axis=1), counts[area].squeeze(), slicer=slicer,
            )
    return (all_samples, missing_from_predictions) if return_missing else all_samples


def _rel(predictions, truth):
    return predictions / (truth + 0.1 * np.std(truth) + 1e-2) - 1


def rmse(predictions, truth):
    return np.sqrt(np.mean(np.power(predictions - truth, 2)))


def std_err(predictions, truth):
    return np.std(predictions - truth, ddof=1)


def mae(predictions, truth):
    return np.mean(np.abs(predictions - truth))


def mad_err(predictions, truth):
    err = predictions - truth
    return np.mean(np.abs(err - mae(predictions, truth)))


def rel_mae(predictions, truth):
    return np.mean(np.abs(_rel(predictions, truth)))


def mad_rel_err(predictions, truth):
    rel = _rel(predictions, truth)
    return np.mean(np.abs(rel - rel_mae(predictions, truth)))


def rel_rmse(predictions, truth):
    return np.sqrt(np.mean(np.power(_rel(predictions, truth), 2)))


def std_rel_err(predictions, truth):
    return np.std(_rel(predictions, truth), ddof=1)


def apply_weekly(metric, sample):
    return [metric(s.predictions, s.truth) for s in sample.values()]


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
                    apply_weekly(metrics_dct[mname], samples),
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
    means = defaultdict(dict)
    for run_name, sample_dct in tqdm(all_samples.items(), desc="Plotting aggregations"):
        weeks = list(next(iter(sample_dct.values())).keys())
        fig, axarr = plt.subplots(len(metrics_dct), 1, constrained_layout=True)
        for ax, (mname, metric) in zip(axarr.flatten(), metrics_dct.items()):
            data = np.row_stack([apply_weekly(metric, s) for s in sample_dct.values()])
            means[metric][run_name] = data.mean(0)
            ax.violinplot(
                np.log(data + 1), showmeans=True, showextrema=False, positions=weeks
            )
            ax.set_xticks(weeks)
            ax.set_ylabel(f"log({mname} + 1)")
            ax.grid(alpha=0.25)
        ax.set_xlabel("Weeks modelled from 1st June, 2020")
        fig.suptitle(f"{run_name}")
        yield fig
        plt.close()

    fig, axarr = plt.subplots(len(metrics_dct), 1, constrained_layout=True)
    for ax, (mname, metric) in zip(axarr.flatten(), metrics_dct.items()):
        for run_name, data in means[metric].items():
            ax.plot(weeks, data, label=run_name)
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

            all_samples[folder][week] = group_samples(
                counts=uk_cases, preds=projections, index="date", by="area"
            )

        all_samples = utils.swaplevel(
            {k: utils.swaplevel(v) for k, v in all_samples.items()}
        )
    if args.pkl_dump is not None:
        with open(args.pkl_dump, "wb") as f:
            pickle.dump(all_samples, f, pickle.HIGHEST_PROTOCOL)

    metrics_dct = {"RMSE": rmse, "MAE": mae}  # "Rel RMSE", "Rel MAE"

    if args.plot_all_areas:
        deck = utils.PdfDeck()
        for fig in plot_all_areas(all_samples, metrics_dct):
            deck.add_figure(fig)
        print("Writing pdf...")
        deck.make(os.path.join(args.outputs_dir, "backtest_evaluation_all_areas.pdf"))

    if args.plot_aggregations:
        all_samples = utils.swaplevel(all_samples)
        deck = utils.PdfDeck()
        for fig in plot_aggregations(all_samples, metrics_dct):
            deck.add_figure(fig)
        deck.figs.insert(0, deck.figs.pop(-1))  # put the aggregate one at the front
        print("Writing pdf...")
        deck.make(
            os.path.join(args.outputs_dir, "backtest_evaluation_aggregations.pdf")
        )
