from covid19_datasets import UKCovid19Data

uk_data = UKCovid19Data()
uk_cases = uk_data.get_cases_data()

# Cases are combined for Cornwall and Isles of Scilly, but they exist as two
# separate Upper Tier Local Authorities:
# E06000052	Cornwall         (population 276,731)
# E06000053	Isles of Scilly  (population 1,089)
# We account both under 'Cornwall' for the purpose of creating a map.

uk_cases.rename(
    index={'Cornwall and Isles of Scilly': 'Cornwall'},
    level='Area name',
    inplace=True)

uk_cases.to_csv('uk_cases.csv')
