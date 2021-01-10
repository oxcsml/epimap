#%%
import pandas as pd
import numpy as np

mapping = {
    "Oxford": "Tehtropolis",
    "Cherwell": "Brynshire",
    "West Oxfordshire": "Bobbingdon",
    "South Oxfordshire": "Hutchintown",
    "Buckinghamshire": "Shehland",
}

order = list(sorted(list(mapping.values())))

areas = pd.read_csv("../../data/areas.csv")
cases = pd.read_csv("../../data/cases.csv")
metadata = pd.read_csv("../../data/metadata.csv")
distances = pd.read_csv("../../data/distances.csv")
nhs_regions = pd.read_csv("../../data/nhs_regions.csv")
traffic_flux = pd.read_csv("../../data/traffic_flux_row-normed.csv")
traffic_flux_transpose = pd.read_csv("../../data/traffic_flux_transpose_row-normed.csv")

areas = areas[areas.area.isin(mapping.keys())]
areas.replace({"area": mapping}, inplace=True)
cases = cases[cases["Area name"].isin(mapping.keys())]
cases.replace({"Area name": mapping}, inplace=True)
cases.drop(columns=["Country"], inplace=True)
cases.rename(columns={"Area name": "area"}, inplace=True)
cases.rename({"AREA": "area"}, inplace=True)
metadata = metadata[metadata.AREA.isin(mapping.keys())]
metadata.replace({"AREA": mapping}, inplace=True)
metadata.rename(columns={"AREA": "area"}, inplace=True)
distances.rename(columns={"Unnamed: 0": "area"}, inplace=True)
traffic_flux.rename(columns={"Unnamed: 0": "area"}, inplace=True)
traffic_flux_transpose.rename(columns={"Unnamed: 0": "area"}, inplace=True)


serial_interval = pd.read_csv("../../data/serial_interval.csv")


def process_df(df):
    df = df[df["area"].isin(mapping.keys())][["area"] + list(mapping.keys())]
    df.replace({"area": mapping}, inplace=True)
    df.rename(columns=mapping, inplace=True)
    df.sort_values(by="area", inplace=True)
    df = df[["area"] + order]

    return df


distances = process_df(distances)
traffic_flux = process_df(traffic_flux)
traffic_flux_transpose = process_df(traffic_flux_transpose)


def norm_flux(df):
    df = df.set_index("area")
    cols = df.columns
    data = df.to_numpy()
    data[np.eye(len(df)).astype(bool)] = 0
    data[np.eye(len(df)).astype(bool)] = 1 - data.sum(axis=1)
    df = pd.DataFrame(data, columns=cols)
    df.index = cols
    return df


traffic_flux = norm_flux(traffic_flux)
traffic_flux_transpose = norm_flux(traffic_flux_transpose)

areas = areas.sort_values(by="area")
cases = cases.sort_values(by="area")
metadata = metadata.sort_values(by="area")

areas.to_csv("data/areas.csv", index=False)
cases.to_csv("data/cases.csv", index=False)
metadata.to_csv("data/metadata.csv", index=False)
distances.to_csv("data/distances.csv", index=False)
nhs_regions.to_csv("data/nhs_regions.csv", index=False)
traffic_flux.to_csv("data/traffic_flux_row-normed.csv")
traffic_flux_transpose.to_csv("data/traffic_flux_transpose_row-normed.csv")
serial_interval.to_csv("data/serial_interval.csv", index=False)

# %%
