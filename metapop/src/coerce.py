from argparse import ArgumentParser

import pandas as pd

if __name__ == "__main__":
    parser = ArgumentParser(
        "Coerce: change output of metapop to inputs like uk_cases.csv"
    )
    parser.add_argument(
        "file", help="path to csv file from metapop to change", type=str
    )
    parser.add_argument("-o", "--output", help="Output file path", type=str)
    args = parser.parse_args()

    start_date = pd.to_datetime("2020-01-30")
    df = pd.read_csv(args.file, header=None).rename(
        index=lambda x: start_date + pd.DateOffset(int(x))
    )
    df.columns.name = "Area name"
    dates = df.index
    df = df.T
    df["Country"] = "nowhere"
    outpath = args.output or "cases.csv"
    df.reset_index(drop=False)[["Country", "Area name", *dates]].to_csv(
            outpath,
            index=False
        )
