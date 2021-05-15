import os
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from matplotlib.dates import date2num as d2n
from collections.abc import Iterable
import numpy as np
month_dict = {"01": "Jan", "02": "Feb", "03": "Mar", "04": "Apr", "05": "May", "06": "Jun", 
            "07": "Jul", "08": "Aug", "09": "Sep", "10": "Oct", "11": "Nov", "12": "Dec"}

def get_xticks_and_xlabels(dfs, col="Date"):
    date_df = pd.concat([df["Date"] for df in dfs])
    dates = date_df.unique()
    months = [month_dict[d.split("-")[1]] for d in dates] 
    date_nums = d2n(dates)
    indx = [d.split("-")[-1] == "01" for d in dates]
    indx = [i for i,x in enumerate(indx) if x]
    date_nums = date_nums[indx]
    months = [months[id] for id in indx]
    return date_nums, months

def subset_df(df, col, vals):
    condition = np.any([df[col]==val for val in vals], axis=0)
    return df.loc[condition]

def make_cases_plot(ax, region, cpred_df, cproj_df, actual_case_df, model_col="royalblue", actual_col="dodgerblue", start="2020-11-30"):
    # ax.set_title(f"{region} cases", fontsize=20)
    # ax.set_ylabel(f"Daily case count", fontsize=15)
    # ax.set_xlabel(f"Date")
    if region == "England":
        england_regions = ['North East and Yorkshire', 'North West', 'Midlands', 'South West', 
                    'East of England', 'South East', 'London']
        reg_actual_df = subset_df(actual_case_df, "NHS_Region", england_regions).groupby("Date", as_index=False).sum()
        reg_cproj_df = subset_df(cproj_df, "area", england_regions).groupby("Date", as_index=False).sum()
        reg_cpred_df = subset_df(cpred_df, "area", england_regions).groupby("Date", as_index=False).sum()
    else:
        reg_actual_df = subset_df(actual_case_df, "NHS_Region", [region])
        reg_cproj_df = subset_df(cproj_df, "area", [region])
        reg_cpred_df = subset_df(cpred_df, "area", [region])
    
    reg_cproj_df = reg_cproj_df[reg_cproj_df["Date"]>start]
    reg_cpred_df = reg_cpred_df[reg_cpred_df["Date"]>start]
    reg_actual_df = reg_actual_df[reg_actual_df["Date"]>start]
    ax.plot(d2n(reg_cpred_df["Date"]), reg_cpred_df["C_50"], color=model_col)
    ax.fill_between(d2n(reg_cpred_df["Date"]), reg_cpred_df["C_25"], reg_cpred_df["C_75"], color=model_col, alpha=0.5)
    ax.fill_between(d2n(reg_cpred_df["Date"]), reg_cpred_df["C_025"], reg_cpred_df["C_975"], color=model_col, alpha=0.25)
    ax.axvline(d2n(reg_cpred_df["Date"])[-1] + 0.5, ls="--", color=model_col, alpha=0.7)

    ax.plot(d2n(reg_cproj_df["Date"]), reg_cproj_df["C_50"], color="k", ls="--")
    ax.fill_between(d2n(reg_cproj_df["Date"]), reg_cproj_df["C_25"], reg_cproj_df["C_75"], color="k", alpha=0.5)
    ax.fill_between(d2n(reg_cproj_df["Date"]), reg_cproj_df["C_025"], reg_cproj_df["C_975"], color="k", alpha=0.25)

    ax.plot(d2n(reg_actual_df["Date"]), reg_actual_df["cases_new"], color=actual_col, alpha=0.4)
    xticks, xlabels = get_xticks_and_xlabels([reg_cpred_df, reg_cproj_df, reg_actual_df])
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=9)

    ax.get_yaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, p: "{:.1f}".format(x/1000)))
    ax.tick_params(axis='y', labelsize=9)

def make_r_plot(ax, region, rt_df, actual_case_df, model_col="royalblue", tick_spacing=0.5, start="2020-11-30"):
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

    reg_rt_df = subset_df(rt_df, "area", [region])
    reg_rt_df = reg_rt_df[reg_rt_df["Date"]>start]
    inferred_df = reg_rt_df.loc[reg_rt_df["provenance"] == "inferred"]
    projected_df = reg_rt_df.loc[reg_rt_df["provenance"] == "projected"]
    ax.axvline(d2n(inferred_df["Date"])[-1] + 0.5, ls="--", color=model_col, alpha=0.7)
    ax.axhline(1.0, ls="--", color="k", alpha=0.9)
    ax.plot(d2n(inferred_df["Date"]), inferred_df["Rt_50"], color=model_col)
    ax.fill_between(d2n(inferred_df["Date"]), inferred_df["Rt_25"], inferred_df["Rt_75"], color=model_col, alpha=0.5)
    ax.fill_between(d2n(inferred_df["Date"]), inferred_df["Rt_025"], inferred_df["Rt_975"], color=model_col, alpha=0.25)

    ax.plot(d2n(projected_df["Date"]), projected_df["Rt_50"], color="k", ls="--")
    ax.fill_between(d2n(projected_df["Date"]), projected_df["Rt_25"], projected_df["Rt_75"], color="k", alpha=0.5)
    ax.fill_between(d2n(projected_df["Date"]), projected_df["Rt_025"], projected_df["Rt_975"], color="k", alpha=0.25)
    actual_case_df = actual_case_df[actual_case_df["Date"] > start]
    xticks, xlabels = get_xticks_and_xlabels([reg_rt_df, actual_case_df])
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=9)
    ax.tick_params(axis='y', labelsize=9)


def create_regional_plot(Rt_file, Cpred_file, Cproj_file, Cactual_file, NHS_regions_file, 
                        save_path="doc/assets/data",
                        regions_to_plot="all"):
    rt_df = pd.read_csv(Rt_file)
    cpred_df = pd.read_csv(Cpred_file)
    cproj_df = pd.read_csv(Cproj_file)
    cactual_df = pd.read_csv(Cactual_file)
    if regions_to_plot=="all":
        regions = ['North East and Yorkshire', 'North West', 'Midlands', 'South West', 
                    'East of England', 'South East', 'London', "England", 'Scotland', 'Wales']
        fig, axs = plt.subplots(4, 5, figsize=(15,8), sharex=True)
        width=5
        height = 4
    elif regions_to_plot=="main_paper":
        regions = ["London", "England", "Scotland", "Wales"]
        fig, axs = plt.subplots(2, 4, figsize=(12,4), sharex=True)
        width=4
        height = 2
    elif regions_to_plot=="appendix":
        regions = ['North East and Yorkshire', 'North West', 'Midlands', 'South West', 
                    'East of England', 'South East']
        fig, axs = plt.subplots(4, 3, figsize=(9,8), sharex=True)
        width=3
        height=4
    else:
        raise ValueError(f"Region selection {regions_to_plot} not found")

    r_id = 0
    plot_id = 0
    for ax_id, ax in enumerate(axs.flatten()):
        region = regions[r_id]
        if plot_id==0:
            make_cases_plot(ax, region, cpred_df, cproj_df, cactual_df)
            if ax_id //width == 0 or ax_id //width == 2:
                ax.set_title(f"{region}", fontsize=14)
            if ax_id % width == 0: 
                ax.set_ylabel(f"Daily cases (K)", fontsize=11)

        if plot_id==1:
            make_r_plot(ax, region, rt_df, cactual_df)
            ax.set_ylim([0., 2.5])

            if ax_id % width == 0:
                ax.set_ylabel(f"Rt", fontsize=11)
        if (ax_id+1) % width == 0:
            plot_id = 1 - plot_id
        r_id = (r_id + 1) % width + width * ((ax_id+1) // (2 * width)) # horrible logic

    for i in range(width):
        ax = axs[height-1, i]
        ax.set_xlabel("Month", fontsize=11)
    fig.tight_layout()
    fig.savefig(save_path + "/" + f"regions_{regions_to_plot}.png")

if __name__ == '__main__':
    args = sys.argv[1:]
    create_regional_plot(*args)
