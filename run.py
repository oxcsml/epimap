import geo.ellipsoid as geo
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import pystan

from datetime import datetime, timedelta
from hashlib import md5
from typing import Dict, Tuple


def plot_areas(case_predictions: pd.DataFrame):
    if not os.path.exists('figs'):
        os.makedirs('figs')

    for area, area_df in case_predictions.groupby(level=0):
        print('Plotting', area)
        day = area_df.index.get_level_values(1)
        lower = area_df['C_lower'].to_numpy()
        median = area_df['C_median'].to_numpy()
        upper = area_df['C_upper'].to_numpy()

        color = [153 / 255.0, 204 / 255.0, 255 / 255.0]
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        ax.fill_between(day, lower, upper, alpha=0.4,
                        color=color, label='2.5% to 97.5% percentiles')
        ax.plot(day, median, lw=4, color=color, label='median')
        ax.legend(loc='upper left', fontsize=12)
        ax.set_title(area, fontsize=12)
        ax.set_xlabel('future days')
        ax.set_ylabel('projected daily cases')
        plt.savefig('figs/' + area + '.pdf')
        plt.close()


def post_process(uk_cases: pd.DataFrame,
                 data: Dict,
                 fit,
                 first_prediction_date: datetime) -> Tuple[pd.DataFrame,
                                                           pd.DataFrame]:
    """Prints model summary statistics and creates CSV file."""

    parameters = ['Ravg', 'gp_length_scale', 'gp_sigma', 'global_sigma',
                  'local_sigma', 'dispersion', 'coupling_rate']
    for p in parameters:
        df = pd.DataFrame(data=fit[p], columns=[p])
        summary = df.describe(percentiles=[0.025, 0.5, 0.975])
        print(summary)
        print()

    df = pd.DataFrame(data=fit['Rt'])
    summary = df.describe(percentiles=[0.025, 0.5, 0.975])
    reproduction_numbers = summary.loc[['2.5%', '50%', '97.5%'], :].transpose()
    reproduction_numbers['area'] = uk_cases.index
    reproduction_numbers = reproduction_numbers.set_index('area')
    reproduction_numbers.columns = ['Rt_lower', 'Rt_median', 'Rt_upper']

    Cproj = fit['Cproj']
    T_proj = data['Tproj']

    cols = ['C_lower', 'C_median', 'C_upper']
    case_predictions = pd.DataFrame(columns=cols + ['area', 'day'])

    for t in range(T_proj):
        data = pd.DataFrame(data=Cproj[:, :, t])
        summary = data.describe(percentiles=[0.025, 0.5, 0.975])
        case_preds = summary.loc[['2.5%', '50%', '97.5%'], :].transpose()
        case_preds.columns = cols
        case_preds['area'] = uk_cases.index
        case_preds['day'] = first_prediction_date + timedelta(days=t)
        case_predictions = case_predictions.append(case_preds)

    case_predictions.set_index(['area', 'day'], inplace=True)

    return reproduction_numbers.sort_index(), case_predictions.sort_index()


def read_data(t_ignore: int=7,
              t_likelihood: int=7,
              t_predict: int=7):
    """Reads UK cases and creates input for Stan model.

    :param t_ignore: The number of most recent days for which the case counts
    are not reliable.
    :param t_likelihood: The number of days for the likelihood to infer Rt.
    :param t_predict: The number of days held out for predictive probabilities
    evaluation.
    """

    infection_profile = pd.read_csv('data/serial_interval.csv')['fit'].to_numpy()
    uk_cases = pd.read_csv('data/uk_cases.csv').set_index('Area name')
    metadata = pd.read_csv('data/metadata.csv').set_index('AREA')

    df = metadata.join(uk_cases, how='inner')
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

    # Use England, Scotland, Wales data.
    df = df[(df['Country'] == 'England') |
            (df['Country'] == 'Scotland') |
            (df['Country'] == 'Wales')]

    non_date_columns = sum([not col.startswith('2020') for col in df.columns])
    # exclude the last 7 days (counts not reliable)
    cases = df.iloc[:, non_date_columns:-t_ignore].to_numpy().astype(int)

    first_prediction_date = datetime.strptime(df.columns[-t_ignore],
                                              '%Y-%m-%d')

    geoloc = df[['LAT', 'LONG']].to_numpy()
    n, _ = geoloc.shape
    geodist = np.zeros((n, n))
    population = df['POPULATION'].to_numpy()

    for i in range(n):
        for j in range(i + 1, n):
            # Vincenty's formula WGS84.
            geodist[i, j] = geo.distance(geoloc[i], geoloc[j]) / 1000
            geodist[j, i] = geodist[i, j]

    # Compute fluxes for radiation model; see "A universal model for mobility
    # and migration patterns" by Simini et al, Nature 484, 96–100(2012).
    flux = np.zeros((n, n))

    for i in range(n):
        source_population = population[i]
        # The distances from the source population, excluding the source.
        sorted_distance_indices = np.argsort(geodist[i, :])[1:]
        # The sorted destination populations.
        sorted_destination_populations = population[sorted_distance_indices]

        # The total population in a circle with radius geodist[i, :], excluding
        # the source population (which is already not in the cumulative sum)
        # and the destination population (which is subtracted).
        population_within_radius = np.cumsum(
            sorted_destination_populations) - sorted_destination_populations

        flx = (source_population * sorted_destination_populations /
               (source_population + population_within_radius) /
               (source_population + sorted_destination_populations +
                population_within_radius))

        flux[i, sorted_distance_indices] = flx
        # According to Simini et al, `flx` should sum to one, but it sums to
        # `1 - source_population / sum(population)`.
        flux[i, i] = 1 - sum(flx)

    t_cond = cases.shape[1] - t_likelihood - t_predict

    data = {
        'N': cases.shape[0],
        'D': len(infection_profile),
        'Tall': cases.shape[1],
        'Tcond': t_cond,
        'Tlik': t_likelihood,
        'Tproj': 21,  # number of days to project forward
        'Count': cases,
        'geoloc': geoloc,
        'geodist': geodist,
        'flux': flux,
        'infprofile': infection_profile
    }

    return df, data, first_prediction_date


def reinflate(reproduction_numbers: pd.DataFrame,
              case_predictions: pd.DataFrame):
    england_map = pd.read_csv('data/england_meta_areas.csv')
    scotland_map = pd.read_csv('data/nhs_scotland_health_boards.csv')
    metadata = pd.read_csv('data/metadata.csv').set_index('AREA')

    # Append a column containing the population ratios of the Upper Tier Local
    # Authorities in the larger regions.
    england_map = england_map.merge(metadata['POPULATION'],
                                    left_on=['Meta area'],
                                    right_on=['AREA'],
                                    how='left') \
                             .merge(metadata['POPULATION'],
                                    left_on=['area'],
                                    right_on=['AREA'],
                                    how='left')
    england_map['ratio'] = (england_map['POPULATION_y'] /
                            england_map['POPULATION_x'])
    england_map = england_map.drop(columns=['POPULATION_x', 'POPULATION_y'])

    scotland_map = scotland_map.merge(metadata['POPULATION'],
                                      left_on=['NHS Scotland Health Board'],
                                      right_on=['AREA'],
                                      how='left') \
                               .merge(metadata['POPULATION'],
                                      left_on=['area'],
                                      right_on=['AREA'],
                                      how='left')
    scotland_map['ratio'] = (scotland_map['POPULATION_y'] /
                             scotland_map['POPULATION_x'])
    scotland_map = scotland_map.drop(columns=['POPULATION_x', 'POPULATION_y'])

    # Reproduction numbers
    england = england_map.merge(reproduction_numbers,
                                left_on=['Meta area'],
                                right_on=['area'],
                                how='left') \
                         .set_index('area') \
                         .drop(columns=['Meta area', 'ratio'])

    scotland = scotland_map.merge(reproduction_numbers,
                                  left_on=['NHS Scotland Health Board'],
                                  right_on=['area'],
                                  how='left') \
                           .set_index('area') \
                           .drop(columns=['NHS Scotland Health Board',
                                          'ratio'])

    reproduction_numbers = reproduction_numbers.append(england) \
                                               .append(scotland)

    # Case predictions
    england = england_map.merge(case_predictions.reset_index(level=1),
                                left_on=['Meta area'],
                                right_on=['area'],
                                how='left')
    england['C_lower'] = england['C_lower'] * england['ratio']
    england['C_median'] = england['C_median'] * england['ratio']
    england['C_upper'] = england['C_upper'] * england['ratio']
    england = england.drop(columns=['Meta area', 'ratio']) \
                     .set_index(['area', 'day'])

    scotland = scotland_map.merge(case_predictions.reset_index(level=1),
                                  left_on=['NHS Scotland Health Board'],
                                  right_on=['area'],
                                  how='left')
    scotland['C_lower'] = scotland['C_lower'] * scotland['ratio']
    scotland['C_median'] = scotland['C_median'] * scotland['ratio']
    scotland['C_upper'] = scotland['C_upper'] * scotland['ratio']
    scotland = scotland.drop(columns=['NHS Scotland Health Board', 'ratio']) \
                       .set_index(['area', 'day'])

    case_predictions = case_predictions.append(england) \
                                       .append(scotland)

    return reproduction_numbers.sort_index(), case_predictions.sort_index()


def stanmodel_cache(model_code, model_name=None):
    """Avoids recompiling the same Stan model if it's already compiled."""
    code_hash = md5(model_code.encode('ascii')).hexdigest()
    if model_name is None:
        cache_fn = 'cached-model-{}.pkl'.format(code_hash)
    else:
        cache_fn = 'cached-{}-{}.pkl'.format(model_name, code_hash)
    try:
        sm = pickle.load(open(cache_fn, 'rb'))
    except:
        sm = pystan.StanModel(model_code=model_code)
        with open(cache_fn, 'wb') as f:
            pickle.dump(sm, f)
    else:
        print('Using cached StanModel')
    return sm, code_hash


def do_modeling(spatialkernel: str='matern12',
                localkernel: str='local',
                globalkernel: str='global',
                metapop: str='radiation_uniform_in',
                observation: str='negative_binomial_3',
                chains: int=4,
                iterations: int=4000,
                pkl_file: str=None):
    """Runs the Stan model, saves samples, and writes predictions to file.

    :param spatialkernel: The kernel to use in the spatial prior GP. Options are
    `matern12`, `matern32`, `matern52`, `exp_quad` or `none`
    :param globalkernel: Whether to include a local kernel. Options are
    `global` or `none`
    :param localkernel: Whether to include a local kernel. Options are
    `local` or `none`
    :param metapop: The metapopulation model for inter-region cross
    infections. Options are `uniform1`, `uniform2` or `none`.
    :param observation: The observation model. Options are `negative_binomial`
    or `poisson`.
    :param chains: The number of MCMC chains.
    :param iterations: The length of each MCMC chain.
    :param pkl_file: The output file for samples.
    """

    # Avoids C++ recompilation if unnecessary in PyStan
    with open('stan_files/Rmap.stan', 'r') as stan_file:
        model_code = stan_file.read()
        model_code = model_code \
            .replace('SPATIAL', spatialkernel) \
            .replace('LOCAL', localkernel) \
            .replace('GLOBAL', globalkernel) \
            .replace('METAPOP', metapop) \
            .replace('OBSERVATION', observation)
    sm, code_hash = stanmodel_cache(model_code)

    # Read the input data.
    uk_cases, data, first_prediction_date = read_data()

    # Use the MD5 hash of the model code if a file name is not provided.
    if pkl_file is None:
        pkl_file = code_hash + '.pkl'

    try:
        # Load samples from the model from a file.
        fit = pickle.load(open(pkl_file, 'rb'))
    except FileNotFoundError:
        # Sample from the model.
        fit = sm.sampling(data=data,
                          iter=iterations,
                          control={'adapt_delta': 0.9},
                          chains=chains)
        with open(pkl_file, 'wb') as f:
            pickle.dump(fit, f, protocol=-1)

    # Post-process the samples.
    reproduction_numbers, case_predictions = post_process(
        uk_cases, data, fit, first_prediction_date)
    reproduction_numbers, case_predictions = reinflate(
        reproduction_numbers, case_predictions)

    # Save the predictions to file.
    if not os.path.exists('projections'):
        os.makedirs('projections')
    reproduction_numbers.to_csv('projections/Rt.csv', float_format='%.5f')
    case_predictions.to_csv('projections/Cproj.csv', float_format='%.5f')

    # Plot the predictions in a PDF for each region.
    plot_areas(case_predictions)


if __name__ == '__main__':
    do_modeling()
