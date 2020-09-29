import pandas as pd
from urllib.request import urlopen
import json

_POPULATION_PATH = 'data/uk_lad_population_estimates.xls'
_TOPO_JSON_URL = 'docs/assets/data/uk_lad_boundaries.json'

population_df = pd.read_excel(_POPULATION_PATH, sheet_name='MYE 5', skiprows=4)
lad_population_df = population_df[['Code', 'Name', 'Estimated Population mid-2019', '2019 people per sq. km', 'Area (sq km)']].rename(columns={
    'Code': 'CODE',
    'Name': 'AREA',
    'Estimated Population mid-2019': 'POPULATION',
    '2019 people per sq. km': 'DENSITY',
    'Area (sq km)': 'SQUARE_KM_AREA'
})

with open(_TOPO_JSON_URL) as f:
    local_authorities_official = json.load(f)

meta_df = pd.DataFrame(
    [(l['properties']['lad20cd'], l['properties']['lad20nm'], l['properties']['lat'], l['properties']['long'])
        for l in local_authorities_official['objects']['Local_Authority_Districts__May_2020__Boundaries_UK_BFC']['geometries']],
    columns=['CODE', 'AREA', 'LAT', 'LONG'])

df = pd.merge(lad_population_df, meta_df, how='right', on=['CODE', 'AREA'])

df = df.set_index('AREA')
df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

# Scotland groups more than one Upper Tier Local Authority into an NHS Health
# Board. We add metadata for the expanded regions.
scotland_map = pd.read_csv('data/nhs_scotland_health_boards.csv').set_index('area')
scotland = scotland_map.join(df, how='inner')
scotland = scotland.groupby('NHS Scotland Health Board').agg(
    {
        'POPULATION': 'sum',
        'SQUARE_KM_AREA': 'sum',
        'LAT': 'mean',
        'LONG': 'mean'
    })
scotland['CODE'] = 'NHS Scotland Health Board'
scotland['DENSITY'] = scotland['POPULATION'] / scotland['SQUARE_KM_AREA']
scotland.index.names = ['AREA']

# England reports very small Upper Tier Local Authorities with their closest
# larger ones. We add metadata for these expanded regions too.
england_map = pd.read_csv('data/england_meta_areas.csv').set_index('area')
england = england_map.join(df, how='inner')
england = england.groupby('Meta area').agg(
    {
        'POPULATION': 'sum',
        'SQUARE_KM_AREA': 'sum',
        'LAT': 'mean',
        'LONG': 'mean'
    })
england['CODE'] = 'England Meta Area'
england['DENSITY'] = england['POPULATION'] / england['SQUARE_KM_AREA']
england.index.names = ['AREA']

df = df.append(scotland).append(england)

df.to_csv('data/metadata.csv')

print('Wrote data/metadata.csv')
