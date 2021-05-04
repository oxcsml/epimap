#%%
import json
import pandas as pd
import numpy as np

regions = [
    "Oxford",
    "Cherwell",
    "South Oxfordshire",
    "Vale of White Horse",
    "West Oxfordshire",
    "Buckinghamshire",
    "West Berkshire",
    "Reading",
    "Swindon",
    "South Northamptonshire",
    "Wokingham",
    "Dacorum",
    "Milton Keynes",
    "Windsor and Maidenhead",
    "Cotswold",
    "Bracknell Forest",
    "Stratford-on-Avon",
    "Slough",
    "Basingstoke and Deane",
    "Three Rivers",
    "Hart",
    "Luton",
    "Watford",
    "Northampton",
    "Central Bedfordshire",
    "Cheltenham",
    "Surrey Heath",
    "Hillingdon",
    "Rushmoor",
    "St Albans",
    "Runnymede",
    "Daventry",
    "Harrow",
    "Warwick",
    "Spelthorne",
    "Wiltshire",
    "Wychavon",
    "Woking",
    "Hertsmere",
    "Hounslow",
    "Tewkesbury",
    "Gloucester",
    "Ealing",
    "Wellingborough",
    "Rugby",
    "Brent",
    "Test Valley",
    "Barnet",
    "Guildford",
    "Welwyn Hatfield",
]

order = sorted(regions)

areas = pd.read_csv("../../data/areas.csv")
cases = pd.read_csv("../../data/cases.csv")
metadata = pd.read_csv("../../data/metadata.csv")
distances = pd.read_csv("../../data/distances.csv")
nhs_regions = pd.read_csv("../../data/nhs_regions.csv")
traffic_flux = pd.read_csv("../../data/traffic_flux_row-normed.csv")
traffic_flux_transpose = pd.read_csv("../../data/traffic_flux_transpose_row-normed.csv")
commute_flow_forward = pd.read_csv("../../data/uk_forward_commute_flow.csv")
commute_flow_reverse = pd.read_csv("../../data/uk_reverse_commute_flow.csv")

areas = areas[areas.area.isin(regions)]
cases = cases[cases["Area name"].isin(regions)]
cases.drop(columns=["Country"], inplace=True)
cases.rename(columns={"Area name": "area"}, inplace=True)
cases.rename({"AREA": "area"}, inplace=True)
metadata = metadata[metadata.AREA.isin(regions)]
metadata.rename(columns={"AREA": "area"}, inplace=True)
distances.rename(columns={"Unnamed: 0": "area"}, inplace=True)
traffic_flux.rename(columns={"Unnamed: 0": "area"}, inplace=True)
traffic_flux_transpose.rename(columns={"Unnamed: 0": "area"}, inplace=True)
commute_flow_forward.rename(columns={"Unnamed: 0": "area"}, inplace=True)
commute_flow_reverse.rename(columns={"Unnamed: 0": "area"}, inplace=True)


serial_interval = pd.read_csv("../../data/serial_interval.csv")


def process_df(df):
    df = df[df["area"].isin(regions)][["area"] + regions]
    df.rename(columns=mapping, inplace=True)
    df.sort_values(by="area", inplace=True)
    df = df[["area"] + order]

    return df


distances = process_df(distances)
traffic_flux = process_df(traffic_flux)
traffic_flux_transpose = process_df(traffic_flux_transpose)
commute_flow_forward = process_df(commute_flow_forward)
commute_flow_reverse = process_df(commute_flow_reverse)


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
commute_flow_forward = norm_flux(commute_flow_forward)
commute_flow_reverse = norm_flux(commute_flow_reverse)

areas = areas.sort_values(by="area")
areas = areas.loc[:, (areas != 0).any(axis=0)]
cases = cases.sort_values(by="area")
metadata = metadata.sort_values(by="area")

areas = areas.loc[:, (areas != 0).any(axis=0)]
remaining_regions = areas.columns[4:]
remaining_regions_ids = [int(name[-1]) for name in remaining_regions]
remaining_regions_index = {
    idx: i + 1
    for i, idx in zip(range(len(remaining_regions_ids)), remaining_regions_ids)
}
nhs_regions = nhs_regions.loc[remaining_regions_ids]

areas.to_csv("tehtropolis/data/areas.csv", index=False)
cases.to_csv("tehtropolis/data/cases.csv", index=False)
metadata.to_csv("tehtropolis/data/metadata.csv", index=False)
distances.to_csv("tehtropolis/data/distances.csv", index=False)
nhs_regions.to_csv("tehtropolis/data/nhs_regions.csv", index=False)
traffic_flux.to_csv("tehtropolis/data/traffic_flux_row-normed.csv")
traffic_flux_transpose.to_csv("tehtropolis/data/traffic_flux_transpose_row-normed.csv")
serial_interval.to_csv("tehtropolis/data/serial_interval.csv", index=False)
commute_flow_forward.to_csv("tehtropolis/data/uk_forward_commute_flow.csv", index=True)
commute_flow_reverse.to_csv("tehtropolis/data/uk_reverse_commute_flow.csv", index=True)

with open("../../data/region-groupings.json", "r") as f:
    groups = json.load(f)

new_groups = {
    remaining_regions_index[int(group)]: [
        mapping[old_area] for old_area in areas if old_area in regions
    ]
    for group, areas in groups.items()
    if int(group) in remaining_regions_ids
}

with open("tehtropolis/data/region-groupings.json", "w") as f:
    json.dump(new_groups, f)
# %%
