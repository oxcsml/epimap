import geo.ellipsoid as geo
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import pystan

from hashlib import md5
from typing import Dict


def plot_areas(out: pd.DataFrame, T_proj: int):

    if not os.path.exists('figs'):
        os.makedirs('figs')

    for index, row in out.iterrows():
        print('Plotting', index)
        lower = row[['C_' + str(t) + '_lower'
                     for t in range(T_proj)]].to_numpy()
        median = row[['C_' + str(t) + '_median'
                      for t in range(T_proj)]].to_numpy()
        upper = row[['C_' + str(t) + '_upper'
                     for t in range(T_proj)]].to_numpy()

        c = [153 / 255.0, 204 / 255.0, 255 / 255.0]
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        ax.fill_between(np.arange(T_proj) + 1, lower, upper, alpha=0.4,
                        color=c, label='2.5% to 97.5% percentiles')
        ax.plot(np.arange(T_proj) + 1, median, lw=4, color=c, label='median')
        ax.legend(loc='upper left', fontsize=12)
        ax.set_title(index, fontsize=12)
        ax.set_xlabel('future days')
        ax.set_ylabel('projected daily cases')
        plt.savefig('figs/' + index + '.pdf')
        plt.close()


def post_process(df: pd.DataFrame,
                 cori_dat: Dict,
                 fit) -> pd.DataFrame:
    """Prints model summary statistics and creates CSV file."""
    parameters = ['Ravg', 'length_scale', 'func_sigma', 'data_sigma',
                  'dispersion', 'immigration_rate']
    for p in parameters:
        data = pd.DataFrame(data=fit[p], columns=[p])
        summary = data.describe(percentiles=[0.025, 0.5, 0.975])
        print(summary)
        print()

    data = pd.DataFrame(data=fit['Rt'])
    summary = data.describe(percentiles=[0.025, 0.5, 0.975])
    reproduction_number = summary.loc[['2.5%', '50%', '97.5%'], :].transpose()
    reproduction_number['area'] = df.index
    reproduction_number = reproduction_number.set_index('area')
    reproduction_number.columns = ['Rtlower', 'Rtmedian', 'Rtupper']

    Cproj = fit['Cproj']
    T_proj = cori_dat['Tproj']

    out = reproduction_number
    for t in range(T_proj):
        data = pd.DataFrame(data=Cproj[:, :, t])
        summary = data.describe(percentiles=[0.025, 0.5, 0.975])
        case_counts = summary.loc[['2.5%', '50%', '97.5%'], :].transpose()
        case_counts['area'] = df.index
        case_counts = case_counts.set_index('area')
        case_counts.columns = ['C_' + str(t) + '_lower',
                               'C_' + str(t) + '_median',
                               'C_' + str(t) + '_upper']
        out = out.join(case_counts, how='inner')

    return out.sort_index()


def read_data():
    """Reads UK cases and creates input for Stan model."""
    infection_profile = pd.read_csv('data/serial_interval.csv')['fit'].to_numpy()
    uk_cases = pd.read_csv('data/uk_cases.csv').set_index('Area name')
    metadata = pd.read_csv('data/metadata.csv').set_index('AREA')

    df = metadata.join(uk_cases, how='inner')
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

    # Use England data.
    df = df[(df['Country'] == 'England') |
            (df['Country'] == 'Scotland') |
            (df['Country'] == 'Wales')]

    non_date_columns = sum([not col.startswith('2020') for col in df.columns])
    # exclude the last 7 days (counts not reliable)
    cases = df.iloc[:, non_date_columns:-7].to_numpy().astype(int)

    geoloc = df[['LAT', 'LONG']].to_numpy()
    n, _ = geoloc.shape
    geodist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            # Vincenty's formula WGS84
            geodist[i, j] = geo.distance(geoloc[i], geoloc[j]) / 1000
            geodist[j, i] = geodist[i, j]

    cori_dat = {
        'N': cases.shape[0],
        'T': cases.shape[1],
        'T0': 7,  # number of days to average over to estimate Rt
        'Tproj': 21,  # number of days to project forward
        'D': len(infection_profile),
        'C': cases,
        'geoloc': geoloc,
        'geodist': geodist,
        'infprofile': infection_profile
    }

    return df, cori_dat


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
    return sm


def do_modeling(kernel: str='exp_quad',
                pkl_file: str='cori-gp-immi-exp_quad_fit.pkl'):
    # Avoids C++ recompilation if unnecessary in PyStan
    with open('cori-gp-immi.stan', 'r') as stan_file:
        model_code = stan_file.read()
        model_code = model_code.replace('KERNEL', kernel)
    sm = stanmodel_cache(model_code)

    # Read the input data.
    df, cori_dat = read_data()

    try:
        fit = pickle.load(open(pkl_file, 'rb'))
    except FileNotFoundError:
        # Sample from the model, or load samples from file.
        fit = sm.sampling(data=cori_dat,
                          iter=4000,
                          control={'adapt_delta': 0.9})
        with open(pkl_file, 'wb') as f:
            pickle.dump(fit, f, protocol=-1)

    # Post-process the samples.
    out = post_process(df, cori_dat, fit)

    out.to_csv('RtCproj-python.csv', float_format='%.5f')
    plot_areas(out, cori_dat['Tproj'])


if __name__ == '__main__':
    do_modeling()
