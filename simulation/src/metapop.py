from argparse import ArgumentParser
import os

import numpy as np
import pandas as pd


def read_csv(pth):
    return pd.read_csv(pth, index_col=0, header=0)


def parse_cmd():
    parser = ArgumentParser(
        "Metapop: deterministic SEIR metapopulation modelling with time varying R"
        " and cross-coupling infection effects."
    )
    parser.add_argument(
        "--a",
        help="Parameter `a` such that mean incubation period is 1/a.",
        required=True,
        type=float,
    )
    parser.add_argument(
        "--gamma",
        help=(
            "Parameter `gamma` such that typical length of"
            " time for which a person is infectious is 1/gamma"
        ),
        type=float,
        required=True,
    )
    parser.add_argument(
        "--init",
        required=True,
        help=("Path to csv containing initial values for S, E, I and R"),
        type=str,
    )
    parser.add_argument(
        "--coupling",
        required=True,
        help=("Path to csv containing coupling constants for cross infection"),
        type=str,
    )
    parser.add_argument(
        "--rt",
        required=True,
        help=(
            "Path to csv containing R_t timeseries by region."
            " The number of rows in this file will determine the"
            " number of time steps for the simulation."
        ),
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help=(
            "Folder in which to save the output files. Defaults to current directory."
        ),
        default=os.getcwd(),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_cmd()
    rt = read_csv(args.rt)
    region_names = rt.columns
    rt = rt.values
    coupling = read_csv(args.coupling).values
    initial = read_csv(args.init).values
    beta = rt * args.gamma

    s = initial[0, :]
    e = initial[1, :]
    i = initial[2, :]
    r = initial[3, :]
    n = s + e + i + r

    susceptible = [s.copy()]
    exposed = [e.copy()]
    infected = [i.copy()]
    recovered = [r.copy()]

    st, et, it, rt = s, e, i, r
    for bt in beta:
        ct = bt * coupling
        sn = st / n
        s = st - sn * np.dot(ct, it)
        e += sn * np.dot(ct, it) - args.a * et
        i += args.a * et - args.gamma * it
        r += args.gamma * it

        st, et, it, rt = s.copy(), e.copy(), i.copy(), r.copy()
        susceptible.append(st)
        exposed.append(et)
        infected.append(it)
        recovered.append(rt)

    filenames = [
        "susceptible.csv",
        "exposed.csv",
        "infected.csv",
        "recovered.csv",
    ]
    all_data = [np.row_stack(k) for k in [susceptible, exposed, infected, recovered]]

    for fname, data in zip(filenames, all_data):
        pd.DataFrame(
                data,
                columns=region_names
                ).to_csv(os.path.join(args.output, fname))
