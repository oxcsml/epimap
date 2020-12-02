import pandas as pd
import json

_POPULATION_PATH = 'data/uk_lad_population_estimates.xls'
_LTLA_TO_NHS_PATH = 'data/ltla_to_nhs.csv'
_TOPO_JSON_URL = 'docs/assets/data/uk_lad_boundaries.json'

population_df = pd.read_excel(_POPULATION_PATH, sheet_name='MYE 5', skiprows=4)
lad_population_df = population_df[['Code', 'Name', 'Estimated Population mid-2019', '2019 people per sq. km', 'Area (sq km)']].rename(columns={
    'Code': 'CODE',
    'Name': 'AREA',
    'Estimated Population mid-2019': 'POPULATION',
    '2019 people per sq. km': 'DENSITY',
    'Area (sq km)': 'SQUARE_KM_AREA'
})

lad_mapping_df = pd.read_csv(_LTLA_TO_NHS_PATH,
                             usecols=['LAU118NM', 'NHS Region']) \
                   .rename(columns={'LAU118NM': 'AREA',
                                    'NHS Region': 'NHS_Region'}) \
                   .set_index('AREA')

with open(_TOPO_JSON_URL) as f:
    local_authorities_official = json.load(f)

meta_df = pd.DataFrame(
    [(l['properties']['lad20cd'], l['properties']['lad20nm'], l['properties']['lat'], l['properties']['long'])
        for l in local_authorities_official['objects']['Local_Authority_Districts__May_2020__Boundaries_UK_BFC']['geometries']],
    columns=['CODE', 'AREA', 'LAT', 'LONG'])

df = pd.merge(lad_population_df, meta_df, how='right', on=['CODE', 'AREA'])

df = df.set_index('AREA')
df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
df = pd.merge(df, lad_mapping_df, how='left', on=['AREA'])

# Scotland groups more than one Upper Tier Local Authority into an NHS Health
# Board. We add metadata for the expanded regions.
scotland_map = pd.read_csv('data/nhs_scotland_health_boards.csv').set_index('area')
scotland = scotland_map.join(df, how='inner')
scotland = scotland.groupby('NHS Scotland Health Board').agg(
    {
        'POPULATION': 'sum',
        'SQUARE_KM_AREA': 'sum',
        'LAT': 'mean',
        'LONG': 'mean',
        'NHS_Region': 'first'
    })
# Corner case: Highlands health board region.
scotland['NHS_Region'] = scotland['NHS_Region'].fillna('Scotland')
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
        'LONG': 'mean',
        'NHS_Region': 'first'
    })
england['CODE'] = 'England Meta Area'
england['DENSITY'] = england['POPULATION'] / england['SQUARE_KM_AREA']
england.index.names = ['AREA']

df = df.append(scotland).append(england)

# Our LTLA to NHS Regions mapping is based on 2018 regions.
# geoportal.statistics.gov.uk has a NUTS-3 to NUTS-2 to NUTS-1 mapping for
# 2018, but I can't find 2020's mapping. We therefore used codes LAU118NM,
# and some regions have changed or were renamed since then. The list below
# sets all the new cases.

df.loc['East Suffolk', 'NHS_Region'] = 'East of England'
df.loc['West Suffolk', 'NHS_Region'] = 'East of England'
df.loc['Buckinghamshire', 'NHS_Region'] = 'South East'
df.loc['Bournemouth, Christchurch and Poole', 'NHS_Region'] = 'South West'
df.loc['Dorset', 'NHS_Region'] = 'South West'
df.loc['Somerset West and Taunton', 'NHS_Region'] = 'South West'
df.loc['Argyll and Bute', 'NHS_Region'] = 'Scotland'
df.loc['Highland', 'NHS_Region'] = 'Scotland'
df.loc['Moray', 'NHS_Region'] = 'Scotland'
df.loc['North Ayrshire', 'NHS_Region'] = 'Scotland'

df.to_csv('data/metadata.csv')

print('Wrote data/metadata.csv')