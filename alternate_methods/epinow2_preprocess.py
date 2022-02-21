from argparse import ArgumentParser
from itertools import starmap
import json

import pandas as pd

if __name__ == "__main__":
    parser = ArgumentParser(description="Change case file format to epiestim style")
    parser.add_argument("in_path", type=str, help="path to input csv")
    parser.add_argument("out_path", type=str, help="where to save outputs")
    parser.add_argument(
        "--region_codes",
        type=str,
        default=None,
        help=(
            "Where to save JSON file with region codes, if at all."
            " This will be dict(area_name -> code)."
        ),
    )
    args = parser.parse_args()

    case_counts = (
        pd.read_csv(args.in_path, index_col=0)
        .set_index("Area name")
        .unstack(0)
        .reset_index(level=1)
        .rename(columns={"Area name": "region", 0: "confirm"})[["confirm", "region"]]
        .rename_axis(index="date")
    )

    case_counts.to_csv(args.out_path)

    if args.region_codes is not None:
        area_codes = dict(
            starmap(lambda k, v: (v, k), enumerate(sorted(case_counts.region.unique())))
        )
        with open(args.region_codes, "w") as f:
            json.dump(area_codes, f)
