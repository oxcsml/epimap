"""Pre-process data for the visualisation website."""

import pandas as pd
import logging
import sys

UK_CASES_PATH = 'data/cases.csv'
UK_METADATA_PATH = 'data/metadata.csv'
OUTPUT_WEBSITE_PATH = 'docs/assets/data/region_site_data.csv'
OUTPUT_DATA_PATH = 'data/region_cases.csv'

# RtCproj_PATH = '../data/RtCproj.csv'

uk_cases = pd.read_csv(UK_CASES_PATH)
metadata = pd.read_csv(UK_METADATA_PATH)

region_df = metadata[["AREA", "NHS_Region"]]
merged_df = uk_cases.merge(region_df, left_on="Area name", right_on="AREA")
merged_df.drop(columns=["AREA", "Area name"], inplace=True)
df = merged_df.groupby(["Country", "NHS_Region"]).sum()
# RtCproj = pd.read_csv(RtCproj_PATH)
df.to_csv(OUTPUT_DATA_PATH, index=True)
print('Wrote', OUTPUT_DATA_PATH)

df = df.stack().to_frame().reset_index().rename(columns={
    'level_2': 'Date',
    0: 'cases_new'
})
df['cases_new_smoothed'] = df['cases_new'].rolling(7, center=True).mean()
df['Date'] = pd.to_datetime(df['Date'])

# Remove the last 5 days of actual cases which are exluded in the modelling due to being unreliable
max_date = df['Date'].max()
df = df[df['Date'] <= max_date - pd.offsets.Day(5)]

df.to_csv(OUTPUT_WEBSITE_PATH, index=False)
print('Wrote', OUTPUT_WEBSITE_PATH)
