import pandas as pd
import logging
import sys
from covid19_datasets import UKCovid19Data

uk_data = UKCovid19Data()
uk_cases = uk_data.get_cases_data()

uk_cases.to_csv('uk_cases.csv')

