import pandas as pd
import numpy as np

# %%

traffic = pd.read_csv("data/uk_traffic.csv", index_col=0)
areas = pd.read_csv("data/areas.csv", index_col=0)

assert (np.diag(traffic) > 0).all()
assert (traffic.sum(axis=1) > 0).all()

# %%

night_populations = areas["population"]
commuting_populations = traffic.sum(axis=1)
noncommuting_populations = night_populations - commuting_populations

# forward flow
forward_flow = traffic + np.diag(noncommuting_populations)
populations = forward_flow.sum(axis=1)
forward_flow = forward_flow.divide(populations, axis="index")
assert np.isclose(forward_flow.sum(axis=1).values, 1).all()
forward_flow.to_csv("data/uk_forward_commute_flow.csv")

# reverse flow
reverse_flow = traffic.T + np.diag(noncommuting_populations)
day_populations = reverse_flow.sum(axis=1)
reverse_flow = reverse_flow.divide(day_populations, axis="index")
assert np.isclose(reverse_flow.sum(axis=1).values, 1).all()
reverse_flow.to_csv("data/uk_reverse_commute_flow.csv")
