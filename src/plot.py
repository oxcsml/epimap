from argparse import ArgumentParser
import os

import matplotlib.pyplot as plt
import pandas as pd


filenames = "susceptible.csv exposed.csv infected.csv recovered.csv".split()


if __name__ == "__main__":
    parser = ArgumentParser("Plot Results")
    parser.add_argument(
        "folder",
        type=str,
        help=(
            "Path to csv file to plot."
            " Must contain files susceptible.csv, exposed.csv, infected.csv and recovered.csv."
            " These files should be of csv type, comma delimited and with the same number of columns (the regions)."
            " The first row will be read as region names."
            " We will assume that there is no index column."
        ),
    )
    args = parser.parse_args()

    outputs = pd.concat(
            {k.replace(".csv", ""): pd.read_csv(os.path.join(args.folder, k))
                for k in filenames
                },
            axis=1
            ).swaplevel(axis=1)

    regions = outputs.columns.levels[0]
    fig, axarr = plt.subplots(len(regions), 1, sharex=True, sharey=True)
    for ax, region in zip(axarr.flat, regions):
        ax.set_title(region, x=0.95, y=0.95, ha="right", va="top")
        outputs[region].plot(ax=ax, legend=False)
        ax.set_ylabel("Population")
        ax.set_xlabel("Timesteps")
        ax.grid(alpha=0.2)
    fig.legend(*ax.get_legend_handles_labels(), ncol=len(regions), bbox_to_anchor=(0.5, 0.05), loc="center")
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)
    plt.show()
