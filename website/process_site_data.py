"""Pre-process data for the visualisation website."""

import pandas as pd
import logging
import sys

UK_CASES_PATH = '../data/uk_cases.csv'
OUTPUT_PATH = 'site_data.csv'
# RtCproj_PATH = '../data/RtCproj.csv'

uk_cases = pd.read_csv(UK_CASES_PATH)
# RtCproj = pd.read_csv(RtCproj_PATH)

df = uk_cases.set_index(['Country', 'Area name']).stack().to_frame().reset_index().rename(columns={
    'Area name': 'area',
    'level_2': 'Date',
    0: 'cases_new'
})
df['cases_new_smoothed'] = df['cases_new'].rolling(7, center=True).mean()
df['Date'] = pd.to_datetime(df['Date'])

# Remove the last 3 days of actual cases which are exluded in the modelling due to being unreliable
max_date = df['Date'].max()
df = df[df['Date'] < max_date - pd.offsets.Day(3)]

df.to_csv(OUTPUT_PATH, index=False)
print('Wrote', OUTPUT_PATH)
