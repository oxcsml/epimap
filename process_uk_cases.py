import numpy as np
from covid19_datasets import UKCovid19Data

uk_data = UKCovid19Data()
uk_cases = uk_data.get_cases_data()

# Scotland has a few negative case counts (usually -1, at most -4, and these
# are very sparse). Set Scotland's negative case counts to 0.
uk_cases.loc[['Scotland'], :] = np.maximum(uk_cases.loc[['Scotland'], :], 0)

# Ignore Scotland's bulk correction on 15 June and take the cases as the mean
# of the day before and after.
uk_cases.loc[['Scotland'], '2020-06-15'] = np.round(
    (uk_cases.loc[['Scotland'], '2020-06-14'] +
     uk_cases.loc[['Scotland'], '2020-06-16']) / 2)

uk_cases.to_csv('data/uk_cases.csv')
