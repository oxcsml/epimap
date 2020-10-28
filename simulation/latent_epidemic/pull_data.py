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
areas = pd.read_csv("../../data/areas.csv")
cases = pd.read_csv("../../data/cases.csv")
metadata = pd.read_csv("../../data/metadata.csv")
distances = pd.read_csv("../../data/distances.csv")
traffic_flux = pd.read_csv("../../data/traffic_flux_row-normed.csv")
traffic_flux_transpose = pd.read_csv("../../data/traffic_flux_transpose_row-normed.csv")

areas = areas[areas.area.isin(mapping.keys())]
areas.replace({"area": mapping}, inplace=True)
cases = cases[cases["Area name"].isin(mapping.keys())]
cases.replace({"Area name": mapping}, inplace=True)
cases.drop(columns=["Country"], inplace=True)
metadata = metadata[metadata.AREA.isin(mapping.keys())]
metadata.replace({"AREA": mapping}, inplace=True)

serial_interval = pd.read_csv("../../data/serial_interval.csv")


def process_df(df):
    df = df[df["Unnamed: 0"].isin(mapping.keys())][
        ["Unnamed: 0"] + list(mapping.keys())
    ]
    df.replace({"Unnamed: 0": mapping}, inplace=True)
    df.rename(mapping, inplace=True)

    return df


distances = process_df(distances)
traffic_flux = process_df(traffic_flux)
traffic_flux_transpose = process_df(traffic_flux_transpose)


def norm_flux(df):
    df = df.set_index("Unnamed: 0")
    cols = df.columns
    data = df.to_numpy()
    data[np.eye(len(df)).astype(bool)] = 0
    data[np.eye(len(df)).astype(bool)] = 1 - data.sum(axis=1)
    df = pd.DataFrame(data, columns=cols)
    df.index = cols
    return df


traffic_flux = norm_flux(traffic_flux)
traffic_flux_transpose = norm_flux(traffic_flux_transpose)

areas.to_csv("data/areas.csv", index=False)
cases.to_csv("data/cases.csv", index=False)
metadata.to_csv("data/metadata.csv", index=False)
distances.to_csv("data/distances.csv", index=False)
traffic_flux.to_csv("data/traffic_flux.csv")
traffic_flux_transpose.to_csv("data/traffic_flux_transpose.csv")
serial_interval.to_csv("data/serial_interval.csv", index=False)
