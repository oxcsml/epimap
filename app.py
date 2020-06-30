import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
from urllib.request import urlopen
import json
import folium

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

GEO_JSON_URL = 'http://geoportal1-ons.opendata.arcgis.com/datasets/687f346f5023410ba86615655ff33ca9_1.geojson'
DATA_PATH = 'Rt.csv'
UK_CASES_PATH = 'uk_cases.csv'

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
with urlopen(GEO_JSON_URL) as response:
  local_authorities_official = json.load(response)

df = pd.read_csv(DATA_PATH)
cases_df = pd.read_csv(UK_CASES_PATH)
cases_df = cases_df.set_index('Area name').drop('Country', axis=1)

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
                    'height': '70vh',
                    'border': 'none'
                }),
    html.H3('UTLA Cases Time Series Plots'),
    dcc.Dropdown(
                id='area-select',
                options=[{'label': a, 'value': a} for a in sorted(cases_df.index)],
                value='Select a UTLA'
            ),
    dcc.Graph(id='area-cases-chart')
])

@app.callback(
    Output('area-cases-chart', 'figure'),
    [Input('area-select', 'value')])
def update_cases_chart(area):
    if area == 'Select a UTLA':
        return

    plot_df = cases_df.loc[area].reset_index().rename(columns={'index': 'Date'})
    plot_df.Date = pd.to_datetime(plot_df.Date)
    plot_df['smoothed'] = plot_df[area].rolling(7, center=True).mean()

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=plot_df.Date, y=plot_df[area], name='New Cases (Actual)',
                            line=dict(color='royalblue', dash='dash')))
    fig.add_trace(go.Scatter(x=plot_df.Date, y=plot_df['smoothed'], name='New Cases (Smoothed)',
                            line=dict(color='royalblue')))
    fig.update_layout(title=f'Daily new cases in {area}',
                    xaxis_title='Date',
                    yaxis_title='Daily New Cases')
    return fig


if __name__ == '__main__':
    app.run_server(debug=True,host = '127.0.0.1')
