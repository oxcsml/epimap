from argparse import ArgumentParser
import os

import numpy as np


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("folder", type=str, help="Where to save the data")
    parser.add_argument("--delimiter", type=str, help="csv delimiter", default=",")
    args = parser.parse_args()

    def save_csv(fname, arr):
        return np.savetxt(
            os.path.join(args.folder, fname), arr, delimiter=args.delimiter, fmt="%d",
        )

    mobility = np.ones((5, 5))
    save_csv("mobility.csv", mobility)

    # initial values
    s = [980, 20, 9, 80, 880]
    e = [10, 880, 990, 900, 20]
    i = [10, 100, 1, 20, 100]
    r = [0, 0, 0, 0, 0]
    init = np.array([s, e, i, r])
    save_csv("initial-values.csv", init)

    # seed epidemic on one place
    s = [400, 1000, 1000, 1000, 1000]
    e = [500, 0, 0, 0, 0]
    i = [100, 0, 0, 0, 0]
    r = [0, 0, 0, 0, 0]
    init_seed = np.array([s, e, i, r])
    save_csv("initial-values-seed.csv", init_seed)
