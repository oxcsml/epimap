import numpy as np


def simulate(
    initial_infections,
    R,
    delay_profile,
    observation_dispersion,
    infectivity_profile,
    infection_dispersion,
    flux_matrix,
    mixing_proportions,
    observation_probability,
    simulation_days,
    infection_distribution="negative_binomial_3",
    observation_model="negative_binomial_3",
):

    N, initial_days = initial_infections.shape
    infective_days = infectivity_profile.shape[0]
    delay_days = delay_profile.shape[0]

    assert (
        initial_days >= infective_days
    ), f"Initial epidemic must be longer than the infectivity profile"

    W_rev = np.flip(infectivity_profile)
    D_rev = np.flip(delay_profile)
    psi_X = infection_dispersion
    psi_C = observation_dispersion
    alpha = mixing_proportions
    F = flux_matrix
    beta = observation_probability

    X = np.zeros((N, simulation_days))
    Z = np.zeros((N, simulation_days))
    X[:, :initial_days] = initial_infections

    for t in range(initial_days, simulation_days):
        Z_t = X[:, t - infective_days : t] @ W_rev
        Z_t = (1 - alpha[t]) * Z_t + alpha[t] * (F @ Z_t)
        Z[:, t] = Z_t
        mu = R[:, t] * Z_t

        if infection_distribution == "negative_binomial_3":
            mu[mu == 0] = 1e-5

            # psi = mu * phi_X
            # var = mu + (mu ** 2 / psi)
            var = (1 + psi_X) * mu
            p = (var - mu) / var
            r = mu ** 2 / (var - mu)

            X[:, t] = np.random.negative_binomial(r, 1 - p)
        elif infection_distribution == "poisson":
            X[:, t] = np.random.poisson(mu)
        else:
            raise ValueError(
                f"{infection_distribution} is not a valid infection distribution"
            )

    C = np.zeros_like(X)
    E = np.zeros_like(X)

    if observation_model == "negative_binomial_3":
        for t in range(initial_days, simulation_days):
            E_t = X[:, t - delay_days : t] @ D_rev
            E[:, t] = E_t
            mu = E_t * beta[t]

            mu[mu == 0] = 1e-5

            # phi = mu * phi_C
            # var = mu + (mu ** 2 / phi)
            var = (1 + psi_C) * mu
            p = (var - mu) / var
            r = mu ** 2 / (var - mu)

            C[:, t] = np.random.negative_binomial(r, 1 - p)
    elif observation_model == "biological":
        raise NotImplementedError()
    else:
        raise ValueError(f"{observation_model} is not a valid observation model")

    return X, C, Z, E