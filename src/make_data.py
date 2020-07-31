from argparse import ArgumentParser
import os

import numpy as np


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("folder", type=str, help="Where to save the data")
    parser.add_argument("--delimiter", type=str, help="csv delimiter", default=",")
    args = parser.parse_args()

    mobility = np.ones((5, 5))
    np.savetxt(
        os.path.join(args.folder, "mobility.csv"), mobility, delimiter=args.delimiter
    )

    # initial values
    s = [980, 20, 9, 80, 880]
    e = [10, 880, 990, 900, 20]
    i = [10, 100, 1, 20, 100]
    r = [0, 0, 0, 0, 0]
    init = np.array([s, e, i, r])
    np.savetxt(
        os.path.join(args.folder, "initial-values.csv"), init, delimiter=args.delimiter
    )
