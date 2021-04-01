import os
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from matplotlib.dates import date2num as d2n

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

def subset_df(df, col, val):
    condition = (df[col] == val)
    return df.loc[condition]

def make_cases_plot(ax, region, cpred_df, cproj_df, actual_case_df, model_col="royalblue", actual_col="dodgerblue"):
    ax.set_title(f"{region} cases")
    ax.set_ylabel(f"Daily case count")
    # ax.set_xlabel(f"Date")
    reg_cpred_df = subset_df(cpred_df, "area", region)
    reg_cproj_df = subset_df(cproj_df, "area", region)
    reg_actual_df = subset_df(actual_case_df, "NHS_Region", region)
    ax.plot(d2n(reg_cpred_df["Date"]), reg_cpred_df["C_50"], color=model_col)
    ax.fill_between(d2n(reg_cpred_df["Date"]), reg_cpred_df["C_25"], reg_cpred_df["C_75"], color=model_col, alpha=0.5)
    ax.fill_between(d2n(reg_cpred_df["Date"]), reg_cpred_df["C_025"], reg_cpred_df["C_975"], color=model_col, alpha=0.25)
    ax.axvline(d2n(reg_cpred_df["Date"])[-1], ls="--", color=model_col, alpha=0.7)

    ax.plot(d2n(reg_cproj_df["Date"]), reg_cproj_df["C_50"], color="k", ls="--")
    ax.fill_between(d2n(reg_cproj_df["Date"]), reg_cproj_df["C_25"], reg_cproj_df["C_75"], color="k", alpha=0.5)
    ax.fill_between(d2n(reg_cproj_df["Date"]), reg_cproj_df["C_025"], reg_cproj_df["C_975"], color="k", alpha=0.25)

    ax.plot(d2n(reg_actual_df["Date"]), reg_actual_df["cases_new"], color=actual_col, alpha=0.4)
    xticks, xlabels = get_xticks_and_xlabels([reg_cpred_df, reg_cproj_df, reg_actual_df])
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)

def make_r_plot(ax, region, rt_df, actual_case_df, model_col="royalblue", tick_spacing=0.5):
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    ax.set_title(f"{region} Rt")
    ax.set_ylabel(f"Rt")
    # ax.set_xlabel(f"Date")
    reg_rt_df = subset_df(rt_df, "area", region)
    
    inferred_df = reg_rt_df.loc[reg_rt_df["provenance"] == "inferred"]
    projected_df = reg_rt_df.loc[reg_rt_df["provenance"] == "projected"]
    ax.axvline(d2n(inferred_df["Date"])[-1], ls="--", color=model_col, alpha=0.7)
    ax.axhline(1.0, ls="--", color="k", alpha=0.9)
    ax.plot(d2n(inferred_df["Date"]), inferred_df["Rt_50"], color=model_col)
    ax.fill_between(d2n(inferred_df["Date"]), inferred_df["Rt_25"], inferred_df["Rt_75"], color=model_col, alpha=0.5)
    ax.fill_between(d2n(inferred_df["Date"]), inferred_df["Rt_2_5"], inferred_df["Rt_97_5"], color=model_col, alpha=0.25)

    ax.plot(d2n(projected_df["Date"]), projected_df["Rt_50"], color="k", ls="--")
    ax.fill_between(d2n(projected_df["Date"]), projected_df["Rt_25"], projected_df["Rt_75"], color="k", alpha=0.5)
    ax.fill_between(d2n(projected_df["Date"]), projected_df["Rt_2_5"], projected_df["Rt_97_5"], color="k", alpha=0.25)

    xticks, xlabels = get_xticks_and_xlabels([reg_rt_df, actual_case_df])
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)

def create_regional_plot(Rt_file, Cpred_file, Cproj_file, Cactual_file, NHS_regions_file, save_dir="doc/assets/data"):
    rt_df = pd.read_csv(Rt_file)
    cpred_df = pd.read_csv(Cpred_file)
    cproj_df = pd.read_csv(Cproj_file)
    cactual_df = pd.read_csv(Cactual_file)
    regions = pd.read_csv(NHS_regions_file, header=0).iloc[:,0].to_list()
    fig, axs = plt.subplots(9, 2, figsize=(15,35))
    ax_list = [axs[r_id, 1] for r_id, _ in enumerate(regions)]
    ax_list[0].get_shared_y_axes().join(*ax_list)
    for r_id, region in enumerate(regions):
        for id in range(2):
            ax = axs[r_id, id]
            if id==0:
                make_cases_plot(ax, region, cpred_df, cproj_df, cactual_df)
            if id==1:
                make_r_plot(ax, region, rt_df, cactual_df)
    fig.savefig(f"{save_dir}/regional_plot.pdf")

if __name__ == '__main__':
    args = sys.argv[1:]
    create_regional_plot(
        Rt_file=args[0],
        Cpred_file=args[1],
        Cproj_file=args[2],
        Cactual_file=args[3],
        NHS_regions_file=args[4],
        save_dir=args[5]
    )
