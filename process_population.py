import pandas as pd

_POPULATION_PATH = 'https://docs.google.com/spreadsheets/d/e/2PACX-1vQ6q_iRClBlQdtE0pFfwWbMa8uTkQVQuTiVqigL09Pv3mEhoUpT0TxsrVC_pxwMhw/pub?gid=106941863&single=true&output=csv'

population_df = pd.read_csv(_POPULATION_PATH, skiprows=6)
population_df = population_df[population_df['AGE GROUP'] == 'All ages'][['CODE', 'AREA', '2018']]
population_df.to_csv('population.csv')

print('Wrote population.csv')