import pandas as pd
from urllib.request import urlopen
import json

_POPULATION_PATH = 'https://raw.githubusercontent.com/rs-delve/Rmap/master/uk_lad_population_estimates.xls?token=ABFV53GAC3YDYVNGGXSXXDS7CBBLW'
_GEO_JSON_URL = 'http://geoportal1-ons.opendata.arcgis.com/datasets/b216b4c8a4e74f6fb692a1785255d777_0.geojson?outSR={%22latestWkid%22:27700,%22wkid%22:27700}'

population_df = pd.read_excel(_POPULATION_PATH, sheet_name='MYE 5', skiprows=4)
utla_population_df = population_df[['Code', 'Name', 'Estimated Population mid-2019', '2019 people per sq. km']].rename(columns={
    'Code': 'CODE',
    'Name': 'AREA',
    'Estimated Population mid-2019': 'POPULATION',
    '2019 people per sq. km': 'DENSITY'
})

# Update Buckinghamshire UTLA to use the county-code instead which is given in the GeoJSON
utla_population_df.CODE = utla_population_df.CODE.replace({
    'E06000060': 'E10000002'
})

with urlopen(_GEO_JSON_URL) as response:
    local_authorities_official = json.load(response)

meta_df = pd.DataFrame(
    [(l['properties']['ctyua19cd'], l['properties']['ctyua19nm'], l['properties']['lat'], l['properties']['long'])
        for l in local_authorities_official['features']],
    columns=['CODE', 'AREA', 'LAT', 'LONG'])

df = pd.merge(utla_population_df, meta_df, how='right', on=['CODE', 'AREA'])
df.to_csv('metadata.csv')

print('Wrote metadata.csv')
