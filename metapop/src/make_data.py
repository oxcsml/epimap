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

    def save_csv(fname, df, fmt, header=True, index=False):
        return df.to_csv(
            os.path.join(args.folder, fname),
            sep=args.delimiter,
            float_format=fmt,
            header=header,
            index=index,
        )

    def save(exp_name):
        def wrapper(func):
            @wraps(func)
            def inner(*args, **kwargs):
                """This enforces saving the cols and index in consistent order"""
                mobility, initial_values, rt = func(*args, **kwargs)
                region_names = initial_values.columns
                save_csv(
                    f"{exp_name}-mobility.csv",
                    mobility[region_names].loc[region_names],
                    fmt="%d",
                    header=False,
                )
                save_csv(
                    f"{exp_name}-initial-values.csv",
                    initial_values[region_names],
                    fmt="%d",
                    header=False,
                )
                save_csv(
                    f"{exp_name}-rt.csv", rt[region_names], fmt="%.5f", header=False
                )
                return

            __register.append(inner)
            return inner

        return wrapper

    @save(exp_name="basic")
    def basic_example():
        region_names = "A B C D E".split()
        mobility = pd.DataFrame(
            np.ones((5, 5)), index=region_names, columns=region_names
        )
        s = [400, 1000, 1000, 1000, 1000]
        e = [500, 0, 0, 0, 0]
        i = [100, 0, 0, 0, 0]
        r = [0, 0, 0, 0, 0]
        initial_values = pd.DataFrame([s, e, i, r], columns=region_names)
        rt = pd.DataFrame(np.full((1000, 5), 2.8), columns=region_names)
        return mobility, initial_values, rt

    @save(exp_name="iid")
    def iid():
        region_names = "A B C D E".split()
        mobility = pd.DataFrame(
            np.zeros((5, 5)), index=region_names, columns=region_names
        )
        s = [800, 800, 800, 800, 800]
        e = [100, 100, 100, 100, 100]
        i = [100, 100, 100, 100, 100]
        r = [0, 0, 0, 0, 0]
        initial_values = pd.DataFrame([s, e, i, r], columns=region_names)
        rt = pd.DataFrame(np.full((10000, 5), 2.8), columns=region_names)
        return mobility, initial_values, rt

    @save(exp_name="identical-uniform-movement")
    def identical_uniform_movement():
        region_names = "A B C D E".split()
        mobility = pd.DataFrame(
            np.ones((5, 5)), index=region_names, columns=region_names
        )
        s = [800, 800, 800, 800, 800]
        e = [100, 100, 100, 100, 100]
        i = [100, 100, 100, 100, 100]
        r = [0, 0, 0, 0, 0]
        initial_values = pd.DataFrame([s, e, i, r], columns=region_names)
        rt = pd.DataFrame(np.full((10000, 5), 2.8), columns=region_names)
        return mobility, initial_values, rt

    @save(exp_name="single")
    def single_region():
        region_names = ["A"]
        mobility = pd.DataFrame(
            np.ones((1, 1)), index=region_names, columns=region_names
        )
        s = [800]
        e = [100]
        i = [100]
        r = [0]
        initial_values = pd.DataFrame([s, e, i, r], columns=region_names)
        rt = pd.DataFrame(np.full((10000, 1), 2.8), columns=region_names)
        return mobility, initial_values, rt

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
        mobility = pd.DataFrame(
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
        return mobility, initial_values, rt

    @save(exp_name="metapop")
    def metapop():
        region_names = [
            "city",
            "town-a",
            "town-b",
            "village-a",
            "village-b",
            "village-c",
        ]
        m = np.array(
            [
                [0, 4000, 3500, 1500, 500, 100],
                [0, 0, 5000, 1000, 100, 10],
                [0, 0, 0, 100, 100, 10],
                [0, 0, 0, 0, 100, 50],
                [0, 0, 0, 0, 0, 50],
                [0, 0, 0, 0, 0, 0],
            ]
        )
        m = np.triu(m)
        mobility = pd.DataFrame(m + m.T, index=region_names, columns=region_names,)
        s = [2e6 - 60, 200000, 150000, 60000, 10000, 1000]
        e = [0, 0, 0, 0, 0, 0]
        i = [60, 0, 0, 0, 0, 0]
        r = [0, 0, 0, 0, 0, 0]
        initial_values = pd.DataFrame([s, e, i, r], columns=region_names)

        city_rt = np.interp(
            np.arange(300), [0, 60, 67, 90, 97, 300], [2.0, 2.0, 0.8, 0.8, 1.15, 1.15],
        )
        town_rt = city_rt * 0.75
        village_rt = city_rt * 0.5
        rts = np.array(
            [
                city_rt,
                town_rt,
                town_rt,
                village_rt,
                village_rt,
                np.full(city_rt.shape, 0.9),
            ]
        ).T
        rt = pd.DataFrame(rts, columns=region_names,)
        return mobility, initial_values, rt

    [f() for f in __register]

    # yw town and village thing
    # put in sensible parameters for stuff from reading paper
