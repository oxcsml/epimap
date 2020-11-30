from argparse import ArgumentParser

import pandas as pd

from metapop import read_csv

"""
get correct case counts
make charts (and save) and put in readme
speak to bobby or someone about running the model on the new data
"""

if __name__ == "__main__":
    parser = ArgumentParser(
        "Coerce: change output of metapop to inputs like uk_cases.csv"
    )
    parser.add_argument(
        "file", help="path to csv file from metapop to change", type=str
    )
    parser.add_argument(
        "--scale",
        help=(
            "constant to scale the case counts with,"
            " for SEIR model you can scale et with a to get case counts"
        ),
        default=1,
        type=float,
    )
    parser.add_argument(
            "--shift",
            help="Shift data along one day",
            action="store_true"
        )
    parser.add_argument("-o", "--output", help="Output file path", type=str)
    args = parser.parse_args()

    start_date = pd.to_datetime("2020-01-30")
    df = (
        read_csv(args.file).rename(
            index=lambda x: start_date + pd.DateOffset(int(x))
        ) * args.scale
    )
    if args.shift:
        df = df.shift(1)
    df.columns.name = "Area name"
    dates = df.index
    df = df.T
    df["Country"] = "UK"
    outpath = args.output or "cases.csv"
    df.reset_index(drop=False)[["Country", "Area name", *dates]].to_csv(
        outpath, index=False
    )
