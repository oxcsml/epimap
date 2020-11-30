#%%
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pandas as pd
from simulation import simulate
from plotting import plot_epidemic
from utils import save_simulation_with_data

append = ""
# append = "simulation/latent_epidemic/"

areas = pd.read_csv(append + "data/areas.csv")
counts = pd.read_csv(append + "data/cases.csv")
metadata = pd.read_csv(append + "data/metadata.csv")
distances = pd.read_csv(append + "data/distances.csv")
radiation_fluxes_01 = pd.read_csv(append + "data/radiation_flux_ls=0.1.csv")
# radiation_fluxes_02 = pd.read_csv(append + "data/radiation_flux_ls=0.2.csv")
# radiation_fluxes_05 = pd.read_csv(append + "data/radiation_flux_ls=0.5.csv")
traffic_flux = pd.read_csv(append + "data/traffic_flux_row-normed.csv")
traffic_flux_transpose = pd.read_csv(
    append + "data/traffic_flux_transpose_row-normed.csv"
)
serial_interval = pd.read_csv(append + "data/serial_interval.csv")

radiation_fluxes_01 = radiation_fluxes_01.to_numpy()[:, 1:].astype(np.float)
# radiation_fluxes_02 = radiation_fluxes_02.to_numpy()[:, 1:].astype(np.float)
# radiation_fluxes_05 = radiation_fluxes_05.to_numpy()[:, 1:].astype(np.float)
traffic_flux = traffic_flux.to_numpy()[:, 1:].astype(np.float)
traffic_flux_transpose = traffic_flux_transpose.to_numpy()[:, 1:].astype(np.float)
serial_interval = serial_interval.to_numpy()[:30, 1].astype(np.float)

N = len(areas)

R_interp = "stepwise"
R_noise = "none"
initial_infection_profile = "one_start"

R_weekly = [
    2.5,
    2.5,
    2.5,
    2.5,
    2.5,
    2.5,
    2.5,
    2.5,
    2.5,
    2.5,
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

if initial_infection_profile == "real":
    initial_infections = counts[counts.columns[1:60]].to_numpy()
elif initial_infection_profile == "one_start":
    initial_infections = counts[counts.columns[1:60]].to_numpy()
    initial_infections[0:2, :] = 0
    initial_infections[3:5, :] = 0
    R_weekly = [2.5, 2.5] + R_weekly
else:
    raise ValueError(
        f"{initial_infection_profile} is not a valid intial infection profile"
    )

days = len(R_weekly) * 7
weeks = np.arange(0, days, 7)
days = np.arange(0, days - 7, 1)

R_weekly = np.array(R_weekly)
R_weekly = np.repeat(R_weekly[np.newaxis, :], N, axis=0)

if R_noise == "multiplicative_pre":
    R_weekly = R_weekly * (1 + 0.1 * np.random.randn(*R_weekly.shape))
elif R_noise == "additive_pre":
    R_weekly = R_weekly + 0.1 * np.random.randn(*R_weekly.shape)

R = np.zeros((N, len(days)))

for i in range(N):
    if R_interp == "linear":
        R[i, :] = interp1d(weeks, R_weekly[i, :], kind="linear")(days)
    elif R_interp == "cubic":
        R[i, :] = interp1d(weeks, R_weekly[i, :], kind="cubic")(days)
    elif R_interp == "stepwise":
        R[i, :] = np.repeat(R_weekly[i, :-1], 7)
    else:
        raise ValueError(f"{R_interp} is not a valid R interpolation")

if R_noise == "multiplicative_post":
    R = R * (1 + 0.1 * np.random.randn(*R.shape))
elif R_noise == "additive_post":
    R = R + 0.1 * np.random.randn(*R.shape)

delay_profile = np.array([0.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.3, 0.1])
N = radiation_fluxes_01.shape[0]
uniform_flux = np.ones((N, N)) * 1 / N
flux_matrices = np.concatenate(
    [
        uniform_flux[:, :, np.newaxis],
        radiation_fluxes_01[:, :, np.newaxis],
        traffic_flux[:, :, np.newaxis],
        traffic_flux_transpose[:, :, np.newaxis],
    ],
    axis=2,
)
# normalise rows
flux_mixing = np.array([0.08, 0.92, 0.0, 0.0])
flux_matrix = flux_matrices @ flux_mixing
mixing_proportions = np.ones_like(R[0, :]) * 0.1  # 0.05 * R[0, :].squeeze() / 2.5
observation_probability = np.ones_like(
    mixing_proportions
)  # observation_probability = np.random.uniform(size=len(mixing_proportions))

observation_dispersion = 0.7
infection_dispersion = 1.5

simulation_days = R.shape[1]

np.random.seed(0)

X, C, Z, E = simulate(
    initial_infections,
    R,
    delay_profile,
    observation_dispersion,
    serial_interval,
    infection_dispersion,
    flux_matrix,
    mixing_proportions,
    observation_probability,
    simulation_days,
    observation_model="biological",
)

# X = X[:, 40:]
plot_epidemic(areas.area, R, X, C, R_interp.capitalize() + " R interpolation")

dates = pd.date_range(counts.columns[1], periods=simulation_days)
area_names = areas.area

X = pd.DataFrame(data=X, index=area_names, columns=dates)
C = pd.DataFrame(data=C, index=area_names, columns=dates)
R = pd.DataFrame(data=R, index=area_names, columns=dates)

params = {
    "delay_profile": list(delay_profile),
    "observation_dispersion": observation_dispersion,
    "infection_dispersion": infection_dispersion,
    "flux_mixing": list(flux_mixing),
    "mixing_proportions": list(mixing_proportions),
}

save_simulation_with_data(X, C, R, params, "test_sim")
# %%
