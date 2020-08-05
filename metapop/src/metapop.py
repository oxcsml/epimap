from argparse import ArgumentParser
import os

import numpy as np

DELIM = ","


def read_csv(pth):
    return np.loadtxt(pth, delimiter=DELIM, dtype=float, ndmin=2)


def parse_cmd():
    parser = ArgumentParser(
        "Metapop: deterministic SEIR metapopulation modelling with time varying R."
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
        "--mobility",
        required=True,
        help=("Path to csv containing mobility patterns"),
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
    mobility = read_csv(args.mobility)
    initial = read_csv(args.init)
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
    population = [n.copy()]

    emigrate = mobility.sum(axis=1)
    mob_colsum = mobility.sum(axis=0)

    st, et, it, rt = s, e, i, r
    for bt in beta:
        sn = s / n
        s = s - bt * it * sn + np.dot(mobility, sn) - emigrate * sn
        e += bt * it * sn - args.a * et + np.dot(mobility, et / n) - emigrate * (et / n)
        i += (
            args.a * et
            - args.gamma * it
            + np.dot(mobility, it / n)
            - emigrate * (it / n)
        )
        r += args.gamma * it + np.dot(mobility, rt / n) - emigrate * (rt / n)
        n += mob_colsum - emigrate

        st, et, it, rt = s.copy(), e.copy(), i.copy(), r.copy()
        susceptible.append(st)
        exposed.append(et)
        infected.append(it)
        recovered.append(rt)
        population.append(n)

    filenames = [
        "susceptible.csv",
        "exposed.csv",
        "infected.csv",
        "recovered.csv",
        "population.csv",
    ]
    all_data = [np.row_stack(k) for k in [susceptible, exposed, infected, recovered, population]]

    for fname, data in zip(filenames, all_data):
        np.savetxt(os.path.join(args.output, fname), data, delimiter=DELIM)
