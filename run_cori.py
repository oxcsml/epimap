import pandas as pd
import pickle
import pystan

from hashlib import md5
from typing import Dict


def post_process(df: pd.DataFrame,
                 cori_dat: Dict,
                 fit):
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

    out = out.sort_index()
    out.to_csv('RtCproj-python.csv', float_format='%.5f')


def read_data():
    """Reads UK cases and creates input for Stan model."""
    infection_profile = pd.read_csv('serial_interval.csv')['fit'].to_numpy()
    uk_cases = pd.read_csv('uk_cases.csv').set_index('Area name')
    metadata = pd.read_csv('metadata.csv').set_index('AREA')

    df = metadata.join(uk_cases, how='inner')
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

    # Use England data.
    df = df[df['Country'] == 'England']
    non_date_columns = sum([not col.startswith('2020') for col in df.columns])

    # exclude the last 7 days (counts not reliable)
    cases = df.iloc[:, non_date_columns:-7].to_numpy().astype(int)

    cori_dat = {
        'N': cases.shape[0],
        'T': cases.shape[1],
        'T0': 7,  # number of days to average over to estimate Rt
        'Tproj': 21,  # number of days to project forward
        'D': len(infection_profile),
        'C': cases,
        'geoloc': df[['LAT', 'LONG']].to_numpy(),
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


def main(resample: bool=True,
         pkl_file: str='cori-gp-immi_fit.pkl'):

    # Avoids C++ recompilation if unnecessary in PyStan
    with open('cori-gp-immi.stan', 'r') as stan_file:
        model_code = stan_file.read()
    sm = stanmodel_cache(model_code)

    # Read the input data.
    df, cori_dat = read_data()

    # Sample from the model, or load samples from file.
    if resample:
        fit = sm.sampling(data=cori_dat,
                          iter=4000,
                          control={'adapt_delta': 0.9})
        with open(pkl_file, 'wb') as f:
            pickle.dump(fit, f, protocol=-1)
    else:
        fit = pickle.load(open(pkl_file, 'rb'))

    post_process(df, cori_dat, fit)


if __name__ == '__main__':
    main(resample=False)
