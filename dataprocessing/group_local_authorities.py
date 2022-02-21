import argparse
import json
from functools import partial
import os

import numpy as np
import pandas as pd

import utils


def read_flux(fpath):
    return pd.read_csv(fpath, index_col=0)


def top_thresh(arr, threshold, lookup=None):
    idx = np.argsort(arr)[::-1]
    lookup = lookup if lookup is not None else arr
    return lookup[idx[: np.argmax((np.cumsum(arr[idx]) >= threshold).values) + 1]]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "flux",
        type=str,
        help=(
            "Proportion of working people who live in area <index> "
            "that work in area <column>."
        ),
    )
    parser.add_argument(
        "flux_t",
        type=str,
        help=(
            "Proportion of people who conduct work in area <index> "
            "that live in area <column>"
        ),
    )
    parser.add_argument("areas", type=str, help="area name to nhs region")
    parser.add_argument("output", type=str, help="json path to save output")
    parser.add_argument("--threshold", type=float, default=0.95)

    args = parser.parse_args()

    flux = read_flux(args.flux)
    flux_t = read_flux(args.flux_t)

    top_exports = flux.apply(
        partial(top_thresh, threshold=args.threshold, lookup=flux.index.values),
        axis=1,
    ).apply(list)

    top_imports = flux_t.apply(
        partial(top_thresh, threshold=args.threshold, lookup=flux.index.values),
        axis=1,
    ).apply(list)

    closest = top_exports.add(top_imports).apply(lambda x: list(dict.fromkeys(x)))

    areas = read_flux(args.areas)
    stub = "nhs_region."
    regions = (
        areas[list(filter(lambda x: x.startswith(stub), areas.columns))]
        .rename(columns=lambda x: int(x.replace(stub, "")))
        .astype(bool)
    )

    def make_groups(regions, closest):
        groups = dict()
        for k in regions.columns:
            in_region = set(regions.index[regions[k]])
            affecting_region = set(closest[in_region].sum(0))
            groups[k] = sorted(in_region.union(affecting_region))
        return groups

    groups = make_groups(regions, closest)

    sizes = pd.Series(utils.map_lowest(len, groups))
    assert sizes.ge(regions.sum(0)).all(), "groups are smaller than the regions"
    print("Group sizes:\n", sizes)

    with open(args.output, "w") as f:
        json.dump(groups, f)
