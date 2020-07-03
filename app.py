import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
from urllib.request import urlopen
import json
import folium

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

GEO_JSON_URL = 'http://geoportal1-ons.opendata.arcgis.com/datasets/687f346f5023410ba86615655ff33ca9_1.geojson'
DATA_PATH = 'Rt.csv'

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
with urlopen(GEO_JSON_URL) as response:
  local_authorities_official = json.load(response)

df = pd.read_csv(DATA_PATH)

for row in df.itertuples():
    match = [l for l in local_authorities_official['features'] if l['properties']['ctyua16nm'] == row.area]
    if match:
        match[0]['properties']['Rt'] = row.Rt
        match[0]['properties']['Area Name'] = match[0]['properties']['ctyua16nm']

m = folium.Map(location=[55.37, 3.4], zoom_start=5)
choropleth = folium.Choropleth(
    geo_data=local_authorities_official,
    data=df,
    columns=['area', 'Rt'],
    key_on='feature.properties.ctyua16nm',
    fill_color='OrRd',
    fill_opacity=0.7,
    line_opacity=0.2,
    legend_name='Rt',
    bins=list((0.0,0.5,0.75,1.0,1.5,2.0,100.0)),
    highlight=True,

)
folium.GeoJsonPopup(fields=['Area Name', 'Rt']).add_to(choropleth.geojson)
choropleth.add_to(m)
map_html = m.save('choropleth.html')

app.layout = html.Div(children=[
    html.H1(children='UK Local Authority Rt Estimates'),
    html.Iframe(srcDoc = open('choropleth.html', 'r').read(), 
                style={
                    'width': '100%',
                    'height': '80vh',
                    'border': 'none'
                })
])

if __name__ == '__main__':
    # app.run_server(debug=True)
    app.run_server(debug=True,host = '127.0.0.1')
