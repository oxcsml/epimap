import pandas as pd

_MOBILITY_PATH = 'https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv'
_OUTPUT_PATH = 'mobility.csv'

df = pd.read_csv(_MOBILITY_PATH)

columns = [
  'sub_region_1', 
  'date', 
  'retail_and_recreation_percent_change_from_baseline',	
  'grocery_and_pharmacy_percent_change_from_baseline',
  'parks_percent_change_from_baseline',	
  'transit_stations_percent_change_from_baseline',	
  'workplaces_percent_change_from_baseline',
  'residential_percent_change_from_baseline']
  
df = df[(df.country_region_code == "GB") & (~df.sub_region_1.isna())][columns].rename(columns={'sub_region_1': 'area'})
df.to_csv('mobility.csv')
