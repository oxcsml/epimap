from argparse import ArgumentParser
import os

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import pandas as pd


filenames = "susceptible.csv exposed.csv infected.csv recovered.csv".split()


def title(ax, region):
    return ax.set_title(region, x=0.95, y=0.9, ha="right", va="top")


def legend(fig, ax):
    lins, labs = ax.get_legend_handles_labels()
    return fig.legend(
        lins, labs, ncol=len(labs), bbox_to_anchor=(0.5, 0.05), loc="center"
    )


def read_csv(fpath, index_col=0):
    return pd.read_csv(fpath, header=0, index_col=index_col)


if __name__ == "__main__":
    parser = ArgumentParser("Plot Results")
    parser.add_argument(
        "folder",
        type=str,
        help=(
            "Path to csv file to plot."
            " Must contain files susceptible.csv, exposed.csv,"
            " infected.csv and recovered.csv."
            " These files should be of csv type, comma delimited"
            " and with the same number of columns (the regions)."
            " The first row will be read as region names and the first column as an index."
        ),
    )
    parser.add_argument(
        "--rt",
        type=str,
        help=(
            "Path to rt csv used for simulation."
            " If given we will plot the R_t timeseries."
        ),
        default=None,
    )
    parser.add_argument(
            "--cases",
            type=str,
            help="Path to case counts file. Expected to be in format of uk_cases.csv",
            default=None
        )
    args = parser.parse_args()

    outputs = pd.concat(
        {
            k.replace(".csv", ""): read_csv(
                os.path.join(args.folder, k)
            )
            for k in filenames
        },
        axis=1,
    ).swaplevel(axis=1)

    regions = outputs.columns.levels[0]

    fig, axarr = plt.subplots(len(regions), 1, sharex=True, sharey=False, squeeze=False)
    for ax, region in zip(axarr.flat, regions):
        title(ax, region)
        outputs[region].plot(ax=ax, legend=False)
        ax.set_ylabel("Population")
        ax.grid(alpha=0.2)
    ax.set_xlabel("Timesteps")
    legend(fig, ax)
    plt.subplots_adjust(hspace=0.05)

    if args.rt is not None:
        rt = read_csv(os.path.join(args.rt))
        npop = outputs.groupby(level=0, axis=1).sum()
        rts = rt * outputs.swaplevel(axis=1)["susceptible"] / npop
        xaxis = outputs.index
        fig, axarr = plt.subplots(len(regions), 1, sharex=True, squeeze=False)
        for ax, region in zip(axarr.flat, regions):
            ax.plot(xaxis, rts[region], label="R_t", zorder=100)
            ax.plot(xaxis[:-1], rt[region], label="R_0", alpha=0.5)
            ax.axhline(1, ls="--", alpha=0.5, label="R_t=1", color="k")
            ax.set_ylabel("Reproduction")
            ax.set_xlabel("Days")
            ax.grid(alpha=0.25)
            title(ax, region)
        legend(fig, ax)
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.1)

    if args.cases is not None:
        cases = read_csv(args.cases, index_col=[0, 1]).T
        cases.index = [k.strftime("%x") for k in pd.to_datetime(cases.index)]
        axarr = cases.plot(subplots=True, xlabel="Date", ylabel="# Cases")
        for ax in axarr.flat:
            ax.grid(alpha=0.25)
        fig = ax.figure
        fig.suptitle("Daily Case Counts", y=0.95)
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.1)

    plt.show()
