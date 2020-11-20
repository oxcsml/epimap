# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"figure.figsize": (15, 12)})

# BASE_DIR = "./fits/backtests/latents_space_1_time_60"
WEEKS = [4, 8, 12, 16]
UK_CASES_DIR = "./data/uk_cases.csv"


class Sample:
    def __init__(self, percentiles, truth, pred_key="c_50", slicer=None):
        self._dates = pd.Index.intersection(percentiles.index, truth.index)
        if slicer is not None:
            self._dates = self._dates[slicer]
        self.percentiles = percentiles.reindex(index=self._dates)
        self.predictions = self.percentiles[pred_key]
        self.truth = truth.reindex(index=self._dates)


def read_csv(pth):
    df = pd.read_csv(pth).rename(columns=str.lower)
    df.date = pd.to_datetime(df.date, format="%Y-%m-%d")
    return df


def group_samples_2(counts, preds, index, by, return_missing=False, slicer=None):
    """Similar to group_samples, but deals with uk_cases format"""
    # grouped_counts = counts.set_index(index).groupby(by)
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
                pred.drop([by], axis=1),
                counts[area].squeeze(),
                slicer=slicer,
            )
    out = (all_samples, missing_from_predictions) if return_missing else all_samples
    return out


### Maybe some of the relative things should be std top over std bottom?


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


def success_codes(predictions, truth):
    """
    Return an integer from {0, 1, 2} depending on how successful the prediction
    was. If true case counts lie within 25% CI of median then return 0,
    if not this but within 47.5% CI of median then return 1, else return 2.
    """
    codify = {"c_025": 2, "c_25": 1, "c_50": 0, "c_75": 0, "c_975": 1}
    return self.percentiles.le(self.truth, axis=0).idxmin(axis=1).map(codify)


def _apply(
    metric, samples
):  # can change metrics later if need to get a different thing from them
    return pd.Series({k: metric(s.predictions, s.truth) for k, s in samples.items()})


def apply(m_dict, samples):
    return pd.concat({k: _apply(m, samples) for k, m in m_dict.items()}, axis=1)


def apply_on_dct(dct, apply_fn, num_levels=1):
    dct_to_return = {}
    for k, v in dct.items():
        if num_levels > 1:
            dct_to_return[k] = apply_on_dct(v, apply_fn, num_levels=num_levels-1)
        else:
            dct_to_return[k] = apply_fn(v)
    return dct_to_return

# load true counts and bring into a nicer format
uk_cases = pd.read_csv(UK_CASES_DIR)
uk_cases = uk_cases.T.drop("Country")
uk_cases.columns = uk_cases.iloc[0]
uk_cases = uk_cases.drop("Area name")
uk_cases.index = pd.to_datetime(uk_cases.index)
uk_cases.index.name = "date"
uk_cases.columns.name = None

# load and format all backtest runs and corresponding true counts

def load_backtest_runs(base_dir):
    backtest_runs = {}
    for weeks in WEEKS:
        proj_path = os.path.join(
            base_dir, f"map_start_2020-06-01_weeks_{weeks}/merged_Cproj.csv"
        )
        proj = read_csv(proj_path)

        # split predictions into first, second and third weeks
        all_samples = {
            1: group_samples_2( 
                counts=uk_cases,
                preds=proj,
                index="date",
                by="area",
                slicer=slice(0, 7),
            ),
            2: group_samples_2(
                counts=uk_cases,
                preds=proj,
                index="date",
                by="area",
                slicer=slice(7, 14),
            ),
            3: group_samples_2(
                counts=uk_cases,
                preds=proj,
                index="date",
                by="area",
                slicer=slice(14, 21),
            ),
            "all": group_samples_2( 
                counts=uk_cases,
                preds=proj,
                index="date",
                by="area",
            ),

        }

        backtest_runs[weeks] = all_samples

    return backtest_runs

from evaluation.areas import area_names

m_dict = {"RMSE": rmse, "Rel RMSE": rel_rmse, "MAE": mae, "Rel MAE": rel_mae}

# %%
apply_fn = lambda samples: apply(m_dict, samples)
metrics = ['RMSE', "MAE"] # "Rel RMSE", "Rel MAE"
BASE_DIRS = [
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

backtest_stats = {
    base_dir: apply_on_dct(load_backtest_runs(os.path.join("./evaluation/Cproj_backtest_space_time/backtests", base_dir)), apply_fn, num_levels=2) for base_dir in BASE_DIRS
}

for area in area_names:
    fig, axes = plt.subplots(len(metrics), 1, figsize=(10, 3 * len(metrics)), constrained_layout=True, sharex=True)
    for metric, ax in zip(metrics, axes):

        for base_dir in BASE_DIRS:
            # backtest_runs = load_backtest_runs(os.path.join("./fits/backtests", base_dir))
            # backtest_area_stats = apply_on_dct(backtest_runs, apply_fn, num_levels=2)
            backtest_area_stats = backtest_stats[base_dir]
            ax.plot(WEEKS, [backtest_area_stats[weeks]['all'][metric][area] for weeks in WEEKS], label=base_dir)
            print(base_dir)
        # ax.plot(WEEKS, [backtest_area_stats[weeks][1][metric][area] for weeks in WEEKS], label="First week")
        # ax.plot(WEEKS, [backtest_area_stats[weeks][2][metric][area] for weeks in WEEKS], label="Second week")
        # ax.plot(WEEKS, [backtest_area_stats[weeks][3][metric][area] for weeks in WEEKS], label="Third week")
        ax.legend(fontsize=7)
        ax.set_ylabel(metric)

    axes[-1].set_xlabel("Weeks modelled from 1st June, 2020")
    fig.suptitle(f"Area: {area}")
    fig.savefig(f"evaluation/plots/model_selection/{area}.pdf")
    plt.close()

    print(area)
