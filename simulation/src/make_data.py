"""
Make coupled metapopulation models with time varying R_t and constant
coupling matrix


In metapop.py we calculate the \phi_ij coupling by r_t * coupling
"""
from argparse import ArgumentParser
from functools import wraps
import os

import numpy as np
import pandas as pd


if __name__ == "__main__":
    __register = list()

    parser = ArgumentParser()
    parser.add_argument("folder", type=str, help="Where to save the data")
    parser.add_argument("--delimiter", type=str, help="csv delimiter", default=",")
    args = parser.parse_args()

    def save_csv(fname, df, header=True, index=True):
        return df.to_csv(
            os.path.join(args.folder, fname),
            sep=args.delimiter,
            header=header,
            index=index,
        )

    def save(exp_name):
        def wrapper(func):
            @wraps(func)
            def inner(*args, **kwargs):
                """This enforces saving the cols and index in consistent order"""
                coupling, initial_values, rt = func(*args, **kwargs)
                region_names = initial_values.columns

                save_csv(
                    f"{exp_name}-coupling.csv",
                    coupling[region_names].loc[region_names],
                )
                save_csv(
                    f"{exp_name}-initial-values.csv",
                    initial_values[region_names],
                )
                save_csv(
                    f"{exp_name}-rt.csv", rt[region_names]
                )
                return

            __register.append(inner)
            return inner

        return wrapper

    @save(exp_name="iid")
    def iid():
        region_names = "A B C D E".split()
        coupling = pd.DataFrame(
            np.eye(5), index=region_names, columns=region_names
        )
        s = [800, 800, 800, 800, 800]
        e = [100, 100, 100, 100, 100]
        i = [100, 100, 100, 100, 100]
        r = [0, 0, 0, 0, 0]
        initial_values = pd.DataFrame([s, e, i, r], columns=region_names)
        rt = pd.DataFrame(np.full((10000, 5), 2.8), columns=region_names)
        coupling = pd.DataFrame(np.eye(5), columns=region_names, index=region_names)
        return coupling, initial_values, rt

    @save(exp_name="identical-uniform-coupling")
    def identical_uniform_coupling():
        region_names = "A B C D E".split()
        coupling = pd.DataFrame(
            np.ones((5, 5)), index=region_names, columns=region_names
        )
        s = [800, 800, 800, 800, 800]
        e = [100, 100, 100, 100, 100]
        i = [100, 100, 100, 100, 100]
        r = [0, 0, 0, 0, 0]
        initial_values = pd.DataFrame([s, e, i, r], columns=region_names)
        rt = pd.DataFrame(np.full((10000, 5), 2.8), columns=region_names)
        coupling = pd.DataFrame(
            np.eye(5) + 0.1*(np.ones(5) - np.eye(5)),
            columns=region_names,
            index=region_names,
        )
        return coupling, initial_values, rt

    @save(exp_name="single")
    def single_region():
        region_names = ["A"]
        coupling = pd.DataFrame(
            np.ones((1, 1)), index=region_names, columns=region_names
        )
        s = [800]
        e = [100]
        i = [100]
        r = [0]
        initial_values = pd.DataFrame([s, e, i, r], columns=region_names)
        rt = pd.DataFrame(np.full((10000, 1), 2.8), columns=region_names)
        return coupling, initial_values, rt

    @save(exp_name="gostic-single")
    def gostic_single_region():
        """
        In the gostic paper it seems that they have some continuous decrease in 
        R, rather than these abrupt changes --- is this due to fewer people
        being susceptible?

        There is a comment in their code simulation.R that describes a variable called 
        CONTINUOUS, but I can't find it being used :s
        """
        region_names = ["A"]
        coupling = pd.DataFrame(
            np.ones((1, 1)), index=region_names, columns=region_names
        )
        s = [2e6 - 60]
        e = [0]
        i = [60]
        r = [0]
        initial_values = pd.DataFrame([s, e, i, r], columns=region_names)
        rt = pd.DataFrame(
            np.interp(
                np.arange(300),
                [0, 60, 67, 90, 97, 300],
                [2.0, 2.0, 0.8, 0.8, 1.15, 1.15],
            ),
            columns=region_names,
        )
        return coupling, initial_values, rt

    @save(exp_name="metapop")
    def metapop():
        region_names = [
            "City",
            "Town",
        ]
        s = [2e6 - 60, 2e5]
        e = [0, 0]
        i = [60, 0]
        r = [0, 0]
        initial_values = pd.DataFrame([s, e, i, r], columns=region_names)

        city_rt = np.interp(
            np.arange(300), [0, 60, 67, 90, 97, 300], [2.0, 2.0, 0.8, 0.8, 1.15, 1.15],
        )
        town_rt = np.interp(
            np.arange(300), [0, 60, 67, 90, 97, 300], [1.0, 1.0, 0.8, 0.8, 1.0, 1.0],
        )
        rts = np.array([city_rt, town_rt,]).T
        rt = pd.DataFrame(rts, columns=region_names)
        coupling = pd.DataFrame(
            np.eye(2) + 0.1*(np.ones(2) - np.eye(2)),
            columns=region_names,
            index=region_names,
        )
        return coupling, initial_values, rt

    [f() for f in __register]
