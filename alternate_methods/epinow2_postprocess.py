from argparse import ArgumentParser
import os
import json

import numpy as np
import pandas as pd
from tqdm import tqdm

RT_PERCENTILES = {
    "Rt_2_5": 2.5,
    "Rt_10": 10,
    "Rt_20": 20,
    "Rt_25": 25,
    "Rt_30": 30,
    "Rt_40": 40,
    "Rt_50": 50,
    "Rt_60": 60,
    "Rt_70": 70,
    "Rt_75": 75,
    "Rt_80": 80,
    "Rt_90": 90,
    "Rt_97_5": 97.5,
}

CASE_PERCENTILES = {
    "C_025": 2.5,
    "C_25": 25,
    "C_50": 50,
    "C_75": 75,
    "C_975": 97.5,
}


def get_area_code(fpath):
    fname = os.path.split(fpath)[-1]
    return int(fname.split("_")[0])


def make_dfs(
    fpaths, region_codes, percs_dct, start_date, weeks_modelled, forecast_days
):

    start_ts = pd.Timestamp(start_date)
    end_ts = (
        start_ts
        + pd.Timedelta(weeks_modelled, unit="W")
        + pd.Timedelta(forecast_days, unit="D")
    )
    dates = pd.date_range(start=start_ts, end=end_ts, freq="D")
    dfs = list()
    for fpath in fpaths:
        code = get_area_code(fpath)
        area_name = region_codes.get(code, code)
        samples = np.loadtxt(fpath)

        days_modelled = samples.shape[1]
        provenance = np.r_[
            np.repeat("inferred", days_modelled - forecast_days),
            np.repeat("projected", forecast_days),
        ]

        df = output_df(
            samples=samples,
            percs_dct=percs_dct,
            dates=dates[-days_modelled:],
            provenance=provenance,
            area=area_name,
        )
        dfs.append(df)
    return pd.concat(dfs, axis=0).sort_index(kind="mergesort")


def output_df(samples, percs_dct, dates, provenance, area):
    """make_df.

    Args:
        samples: Samples
        percs_dct:
        dates:
        provenance:
        area:
    """
    percentiles = np.percentile(
        samples, list(percs_dct.values()), interpolation="linear", axis=0
    ).transpose()

    df = pd.DataFrame(percentiles, columns=list(percs_dct.keys()))
    df.insert(0, "Date", dates)
    df["provenance"] = provenance
    df.index = pd.Index([area] * len(dates), name="area")
    return df


if __name__ == "__main__":
    parser = ArgumentParser("Postprocess epinow2 samples")

    parser.add_argument(
        "input_folder",
        type=str,
        help=(
            "Path to folder from which to read .txt files of samples. "
            "Filenames should be <area code/name>_<r/case>_samples.txt. "
            "Columns should be times and rows samples."
        ),
    )
    parser.add_argument(
        "output_folder", type=str, help="Folder in which to save the output csvs"
    )
    parser.add_argument("first_day_modelled", type=str, help="To begin date stamps")
    parser.add_argument("weeks_modelled", type=int, help="Number of weeks modelled")
    parser.add_argument("forecast_days", type=int, help="Number of days of forecast")
    parser.add_argument(
        "--prefix", type=str, help="Filename prefix for output files", default=""
    )
    parser.add_argument(
        "--region_codes",
        type=str,
        default=None,
        help=("Path to JSON file with region codes." "Format area_name -> code."),
    )
    parser.add_argument("--progress", action="store_true", help="Show progress")

    args = parser.parse_args()

    if args.region_codes is not None:
        with open(args.region_codes, "r") as f:
            region_codes = json.loads(f.read())
            region_codes = {v: k for k, v in region_codes.items()}
    else:
        region_codes = {}

    def get_fpaths(suffix, filenames):
        return map(
            lambda x: os.path.join(args.input_folder, x),
            filter(lambda x: x.endswith(suffix), filenames),
        )

    filenames = sorted(os.listdir(args.input_folder))
    r_fpaths = get_fpaths("r_samples.txt", filenames)
    case_fpaths = get_fpaths("case_samples.txt", filenames)

    if args.progress:
        r_fpaths = tqdm(list(r_fpaths), desc="Processing Rt")
    rt = make_dfs(
        fpaths=r_fpaths,
        region_codes=region_codes,
        percs_dct=RT_PERCENTILES,
        start_date=args.first_day_modelled,
        weeks_modelled=args.weeks_modelled,
        forecast_days=args.forecast_days
    )
    rt.to_csv(os.path.join(args.output_folder, f"{args.prefix}_Rt.csv"))

    if args.progress:
        case_fpaths = tqdm(list(case_fpaths), desc="Processing Cases")
    cases_df = make_dfs(
        fpaths=case_fpaths,
        region_codes=region_codes,
        percs_dct=CASE_PERCENTILES,
        start_date=args.first_day_modelled,
        weeks_modelled=args.weeks_modelled,
        forecast_days=args.forecast_days
    )
    cases_df.query("provenance=='inferred'").to_csv(
        os.path.join(args.output_folder, f"{args.prefix}_Cpred.csv")
    )
    cases_df.query("provenance=='projected'").to_csv(
        os.path.join(args.output_folder, f"{args.prefix}_Cproj.csv")
    )
