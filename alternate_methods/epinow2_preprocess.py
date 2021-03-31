from argparse import ArgumentParser

import pandas as pd

if __name__ == "__main__":
    parser = ArgumentParser(description="Change case file format to epiestim style")
    parser.add_argument("in_path", type=str, help="path to input csv")
    parser.add_argument("out_path", type=str, help="where to save outputs")
    args = parser.parse_args()

    (
        pd.read_csv(args.in_path, index_col=0)
        .set_index("Area name")
        .unstack(0)
        .reset_index(level=1)
        .rename(columns={"Area name": "region", 0: "confirm"})[["confirm", "region"]]
        .rename_axis(index="date")
        .to_csv(args.out_path)
    )
