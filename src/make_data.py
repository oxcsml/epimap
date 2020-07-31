from argparse import ArgumentParser
import os

import numpy as np


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("folder", type=str, help="Where to save the data")
    args = parser.parse_args()

    mobility = np.ones((5, 5))
    np.savetxt(os.path.join(args.folder, "mobility.csv"), mobility, delimiter=',')

