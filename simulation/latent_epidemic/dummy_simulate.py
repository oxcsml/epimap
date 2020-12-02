#%%
import numpy as np
from simulation import simulate
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# %%

N = 10

R = [
    2.5,
    2.5,
    2.0,
    0.7,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1.0,
    1.1,
    1.2,
    1.3,
    1.3,
    1.3,
    1.4,
    1.4,
    1.4,
]
days = len(R) * 7
weeks = np.arange(0, days, 7)
days = np.arange(0, days - 7, 1)

R_interp = "linear"

if R_interp == "linear":
    R = interp1d(weeks, R, kind="linear")(days)
elif R_interp == "cubic":
    R = interp1d(weeks, R, kind="cubic")(days)
elif R_interp == "stepwise":
    R = np.repeat(R[:-1], 7)

R = np.repeat(R[np.newaxis, :], N, axis=0)

initial_infections = 10 * np.repeat(
    np.array([1, 2, 4, 6, 10, 15])[np.newaxis, :], N, axis=0
)
delay_profile = np.array([0.1, 0.2, 0.3, 0.2, 0.1, 0.1])
observation_dispersion = 250
infectivity_profile = np.array([0.1, 0.2, 0.3, 0.2, 0.1, 0.1])
infection_dispersion = 30
flux_matrix = np.ones((N, N)) * 1 / N
mixing_proportions = np.ones(len(days)) * 0.05
observation_probability = np.ones(len(days)) * 1.0
simulation_days = len(days)

epi_scale = 1.0

X, C = simulate(
    initial_infections / epi_scale,
    R,
    delay_profile,
    observation_dispersion,
    infectivity_profile,
    infection_dispersion,
    flux_matrix,
    mixing_proportions,
    observation_probability * epi_scale,
    simulation_days,
    # infection_distribution="poisson",
)

# X = X[:, 40:]

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

ax1.plot(X.T, color="tab:blue", alpha=0.3)
ax1.plot(C.T, color="tab:orange", alpha=0.3)
ax1.set_ylabel("Case numbers")

ax2.plot(R[0, :], color="tab:red")
ax2.set_xlabel("Days into epidemic")
ax2.set_ylabel("R")

import matplotlib.patches as mpatches

handles = [
    mpatches.Patch(color=c, label=k)
    for (k, c) in {
        "Latent Epidemic": "tab:blue",
        "Observed Cases": "tab:orange",
    }.items()
]
ax1.legend(handles=handles)

plt.suptitle(R_interp.capitalize() + " R interpolation")