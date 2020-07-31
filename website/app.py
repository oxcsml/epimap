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
import numpy as np

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

GEO_JSON_URL = 'http://geoportal1-ons.opendata.arcgis.com/datasets/687f346f5023410ba86615655ff33ca9_1.geojson'
DATA_PATH = 'fits/Rmap-exp_quad-none-negative_binomial_RtCproj.csv'
print(DATA_PATH)
UK_CASES_PATH = 'data/uk_cases.csv'

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
with urlopen(GEO_JSON_URL) as response:
  local_authorities_official = json.load(response)

df = pd.read_csv(DATA_PATH)
cases_df = pd.read_csv(UK_CASES_PATH)
cases_df = cases_df.set_index('Area name').drop('Country', axis=1)


minRt = 0.5 # np.floor(min(df.Rtmedian)*2.0)/2.0
maxRt = 1.5 # np.ceil(max(df.Rtmedian)*2.0)/2.0
Rtbins = np.exp(np.linspace(np.log(minRt),np.log(maxRt),9))
print(Rtbins)

df.Rtmedian = np.minimum(maxRt,np.maximum(minRt,df.Rtmedian))

minCproj = np.ceil(min(df.Cprojmedian))
maxCproj = np.ceil(max(df.Cprojmedian))
Cprojbins = np.concatenate((np.array([0.0]), np.exp(np.linspace(np.log(minCproj),np.log(maxCproj),6))))
print(Cprojbins)

for row in df.itertuples():
    match = [l for l in local_authorities_official['features'] if l['properties']['ctyua16nm'] == row.area]
    if match:
        match[0]['properties']['Rtlower'] = row.Rtlower
        match[0]['properties']['Rtmedian'] = row.Rtmedian
        match[0]['properties']['Rtupper'] = row.Rtupper
        match[0]['properties']['Cprojlower'] = row.Cprojlower
        match[0]['properties']['Cprojmedian'] = row.Cprojmedian
        match[0]['properties']['Cprojupper'] = row.Cprojupper
        match[0]['properties']['Area Name'] = match[0]['properties']['ctyua16nm']

Rtmap = folium.Map(location=[52.9, -2.0], zoom_start=7)
Rtchoropleth = folium.Choropleth(
    geo_data=local_authorities_official,
    data=df,
    columns=['area', 'Rtmedian'],
    key_on='feature.properties.ctyua16nm',
    fill_color='OrRd',
    fill_opacity=0.7,
    line_opacity=0.2,
    legend_name='Rtmedian',
    bins=list(Rtbins),
    # bins=list((0.0,0.5,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.5)),
    highlight=True,

)
folium.GeoJsonPopup(fields=['Area Name','Rtlower','Rtmedian','Rtupper']).add_to(Rtchoropleth.geojson)
Rtchoropleth.add_to(Rtmap)
Rtmap_html = Rtmap.save('Rtmap.html')

Cprojmap = folium.Map(location=[52.9, -2.0], zoom_start=7)
Cprojchoropleth = folium.Choropleth(
    geo_data=local_authorities_official,
    data=df,
    columns=['area', 'Cprojmedian'],
    key_on='feature.properties.ctyua16nm',
    fill_color='OrRd',
    fill_opacity=0.7,
    line_opacity=0.2,
    legend_name='Cprojmedian',
    # bins=list((0.0,1.0,2.0,5.0,10.0,20.0,50.0)),
    bins=list(Cprojbins),
    highlight=True,

)
folium.GeoJsonPopup(fields=['Area Name','Cprojlower','Cprojmedian','Cprojupper']).add_to(Cprojchoropleth.geojson)
Cprojchoropleth.add_to(Cprojmap)
Cprojmap_html = Cprojmap.save('Cprojmap.html')

app.layout = html.Div(children=[
    html.H1(children=DATA_PATH),
    html.H2(children='UKLA Rt and projected case counts estimates'),
    html.Iframe(srcDoc = open('Rtmap.html', 'r').read(), 
                style={
                    'width': '49%',
                    'height': '600pt',
                }),
    html.Iframe(srcDoc = open('Cprojmap.html', 'r').read(), 
                style={
                    'width': '49%',
                    'height': '600pt',
                }),
    html.H3('UTLA Cases Time Series Plots'),
    dcc.Dropdown(
                id='area-select',
                options=[{'label': a, 'value': a} for a in sorted(cases_df.index)],
                value='Lewisham'
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
    # plot_df['smoothed'] = plot_df[area].rolling(7, center=True).mean()

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=plot_df.Date, y=plot_df[area], name='New Cases (Actual)',
                            line=dict(color='royalblue', dash='dash')))
    # fig.add_trace(go.Scatter(x=plot_df.Date, y=plot_df['smoothed'], name='New Cases (Smoothed)',
    #                         line=dict(color='royalblue')))
    fig.update_layout(title=f'Daily new cases in {area}',
                    xaxis_title='Date',
                    yaxis_title='Daily New Cases')
    return fig


if __name__ == '__main__':
    # app.run_server(debug=True, port=8050)
    app.run_server(debug=True,host = '127.0.0.1')
