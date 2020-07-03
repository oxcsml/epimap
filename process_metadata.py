import pandas as pd
from urllib.request import urlopen
import json

_POPULATION_PATH = 'https://docs.google.com/spreadsheets/d/e/2PACX-1vQ6q_iRClBlQdtE0pFfwWbMa8uTkQVQuTiVqigL09Pv3mEhoUpT0TxsrVC_pxwMhw/pub?gid=106941863&single=true&output=csv'
_GEO_JSON_URL = 'http://geoportal1-ons.opendata.arcgis.com/datasets/687f346f5023410ba86615655ff33ca9_1.geojson'

population_df = pd.read_csv(_POPULATION_PATH, skiprows=6)
population_df = population_df[population_df['AGE GROUP']
                              == 'All ages'][['CODE', 'AREA', '2018']]

with urlopen(_GEO_JSON_URL) as response:
    local_authorities_official = json.load(response)

lat_long_df = pd.DataFrame(
    [(l['properties']['ctyua16cd'], l['properties']['ctyua16nm'], l['properties']['lat'], l['properties']['long'])
        for l in local_authorities_official['features']],
    columns=['CODE', 'AREA', 'LAT', 'LONG'])

df = pd.merge(population_df, lat_long_df, how='right', on=['CODE', 'AREA']).rename(columns={'2018': 'POPULATION'})
df.to_csv('metadata.csv')

print('Wrote metadata.csv')
