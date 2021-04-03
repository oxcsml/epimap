import pandas as pd
import sys


def process_csvs(input_r_filename: str='fits/0_Rt.csv',
                 input_exceed_filename: str='fits/0_Pexceed.csv',
                 input_weekly_filename: str='fits/0_Cweekly.csv',
                 input_Cpred_filename: str='fits/0_Cpred.csv',
                 input_Cproj_filename: str='fits/0_Cproj.csv',
                 input_Xpred_filename: str='fits/0_Xpred.csv',
                 input_Xproj_filename: str='fits/0_Xproj.csv',
                 output_r_filename: str='website/default/Rt.csv',
                 output_exceed_filename: str='website/default/Pexceed.csv',
                 output_weekly_filename: str='website/default/Cweekly.csv',
                 output_Cpred_filename: str='website/default/Cpred.csv',
                 output_Cproj_filename: str='website/default/Cproj.csv',
                 output_Xpred_filename: str='website/default/Xpred.csv',
                 output_Xproj_filename: str='website/default/Xproj.csv'
):

    reproduction_numbers = pd.read_csv(input_r_filename) \
                             .rename(columns={'Date': 'day'})
    reproduction_numbers['day'] = pd.to_datetime(reproduction_numbers['day'],
                                                 format='%Y-%m-%d')
    reproduction_numbers = reproduction_numbers.set_index(['area', 'day'])

    exceedance_probs = pd.read_csv(input_exceed_filename) \
                         .rename(columns={'Date': 'day'})
    exceedance_probs['day'] = pd.to_datetime(exceedance_probs['day'],
                                             format='%Y-%m-%d')
    exceedance_probs = exceedance_probs.set_index(['area', 'day'])

    weekly = pd.read_csv(input_weekly_filename) \
                    .rename(columns={'Date': 'day'})
    weekly['day'] = pd.to_datetime(weekly['day'], format='%Y-%m-%d')
    weekly = weekly.set_index(['area', 'day'])


    Cpredictions = pd.read_csv(input_Cpred_filename) \
                    .rename(columns={'Date': 'day'})
    Cpredictions['day'] = pd.to_datetime(Cpredictions['day'], format='%Y-%m-%d')
    Cpredictions = Cpredictions.set_index(['area', 'day'])

    Cprojections = pd.read_csv(input_Cproj_filename) \
                    .rename(columns={'Date': 'day'})
    Cprojections['day'] = pd.to_datetime(Cprojections['day'], format='%Y-%m-%d')
    Cprojections = Cprojections.set_index(['area', 'day'])

    Xpredictions = pd.read_csv(input_Xpred_filename) \
                    .rename(columns={'Date': 'day'})
    Xpredictions['day'] = pd.to_datetime(Xpredictions['day'], format='%Y-%m-%d')
    Xpredictions = Xpredictions.set_index(['area', 'day'])

    Xprojections = pd.read_csv(input_Xproj_filename) \
                    .rename(columns={'Date': 'day'})
    Xprojections['day'] = pd.to_datetime(Xprojections['day'], format='%Y-%m-%d')
    Xprojections = Xprojections.set_index(['area', 'day'])


    (reproduction_numbers, exceedance_probs, weekly, 
     Cpredictions, Cprojections, Xpredictions, Xprojections) = reinflate(
     reproduction_numbers, exceedance_probs, weekly, 
     Cpredictions, Cprojections, Xpredictions, Xprojections)

    reproduction_numbers.index.names = ['area', 'Date']
    exceedance_probs.index.names = ['area', 'Date']
    weekly.index.names = ['area', 'Date']
    Cpredictions.index.names = ['area', 'Date']
    Cprojections.index.names = ['area', 'Date']
    Xpredictions.index.names = ['area', 'Date']
    Xprojections.index.names = ['area', 'Date']

    reproduction_numbers.to_csv(output_r_filename, float_format='%.2f')
    exceedance_probs.to_csv(output_exceed_filename, float_format='%.2f')
    weekly.to_csv(output_weekly_filename, float_format='%3d')
    Cpredictions.to_csv(output_Cpred_filename, float_format='%2.1f')
    Cprojections.to_csv(output_Cproj_filename, float_format='%2.1f')
    Xpredictions.to_csv(output_Xpred_filename, float_format='%2.1f')
    Xprojections.to_csv(output_Xproj_filename, float_format='%2.1f')


def reinflate(reproduction_numbers: pd.DataFrame,
              exceedance_probs: pd.DataFrame,
              weekly: pd.DataFrame,
              Cpredictions: pd.DataFrame,
              Cprojections: pd.DataFrame,
              Xpredictions: pd.DataFrame,
              Xprojections: pd.DataFrame,
              divide_predictions_by_population_ratio: bool=False):

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

    # `Highland` appears as both an UTLA area and an NHS Scotland Health Board
    # name, and some extra care is required to get the population ratios
    # correct for the areas in the `Highland` health board.
    scotland_map = scotland_map.merge(metadata[['POPULATION', 'CODE']],
                                      left_on=['NHS Scotland Health Board'],
                                      right_on=['AREA'],
                                      how='left')
    scotland_map = scotland_map[scotland_map.CODE ==
                                'NHS Scotland Health Board'] \
                   .drop(columns=['CODE'])

    scotland_map = scotland_map.merge(metadata[metadata.CODE !=
                                               'NHS Scotland Health Board']
                                              ['POPULATION'],
                                      left_on=['area'],
                                      right_on=['AREA'],
                                      how='left')

    scotland_map['ratio'] = (scotland_map['POPULATION_y'] /
                             scotland_map['POPULATION_x'])
    scotland_map = scotland_map.drop(columns=['POPULATION_x', 'POPULATION_y'])

    # Reproduction numbers
    # `Highland` is both the name of an NHS Scotland Health Board and an area
    # in the health board. To avoid a duplicate area key, we drop it as a
    # health board, and keep it as a local area (in this case the two rows
    # would have been equal in any event).
    reproduction_numbers = reproduction_numbers.rename(
        index={'Highland': 'Highland (NHS Scotland Health Board)'})


    england = england_map.merge(reproduction_numbers.reset_index(level=1),
                                left_on=['Meta area'],
                                right_on=['area'],
                                how='left') \
                         .set_index(['area', 'day']) \
                         .drop(columns=['Meta area', 'ratio'])

    scotland = scotland_map.merge(reproduction_numbers.reset_index(level=1),
                                  left_on=['NHS Scotland Health Board'],
                                  right_on=['area'],
                                  how='left') \
                           .set_index(['area', 'day']) \
                           .drop(columns=['NHS Scotland Health Board',
                                          'ratio'])

    reproduction_numbers = reproduction_numbers.append(england) \
                                               .append(scotland)

    # Exceedance probabilities
    # `Highland` is both the name of an NHS Scotland Health Board and an area
    # in the health board. To avoid a duplicate area key, we drop it as a
    # health board, and keep it as a local area (in this case the two rows
    # would have been equal in any event).
    exceedance_probs = exceedance_probs.rename(
        index={'Highland': 'Highland (NHS Scotland Health Board)'})


    england = england_map.merge(exceedance_probs.reset_index(level=1),
                                left_on=['Meta area'],
                                right_on=['area'],
                                how='left') \
                         .set_index(['area', 'day']) \
                         .drop(columns=['Meta area', 'ratio'])

    scotland = scotland_map.merge(exceedance_probs.reset_index(level=1),
                                  left_on=['NHS Scotland Health Board'],
                                  right_on=['area'],
                                  how='left') \
                           .set_index(['area', 'day']) \
                           .drop(columns=['NHS Scotland Health Board',
                                          'ratio'])

    exceedance_probs = exceedance_probs.append(england) \
                                       .append(scotland)
    # Weekly Cases 
    england = england_map.merge(weekly.reset_index(level=1),
                                left_on=['Meta area'],
                                right_on=['area'],
                                how='left')

    if divide_predictions_by_population_ratio:
        england['C_weekly'] = england['C_weekly'] * england['ratio']

    england = england.drop(columns=['Meta area', 'ratio']) \
                     .set_index(['area', 'day'])

    scotland = scotland_map.merge(weekly.reset_index(level=1),
                                  left_on=['NHS Scotland Health Board'],
                                  right_on=['area'],
                                  how='left')

    if divide_predictions_by_population_ratio:
        scotland['C_weekly'] = scotland['C_weekly'] * scotland['ratio']

    scotland = scotland.drop(columns=['NHS Scotland Health Board', 'ratio']) \
                       .set_index(['area', 'day'])

    # We only keep the case for `Highland` as a local area, not
    # weekly = weekly.drop(['Highland'])
    if divide_predictions_by_population_ratio:
        weekly = weekly.rename(
            index={'Highland': 'Highland (NHS Scotland Health Board)'})

    weekly = weekly.append(england) \
                   .append(scotland)

    # Case predictions 
    england = england_map.merge(Cpredictions.reset_index(level=1),
                                left_on=['Meta area'],
                                right_on=['area'],
                                how='left')

    if divide_predictions_by_population_ratio:
        england['C_025'] = england['C_025'] * england['ratio']
        england['C_10'] = england['C_10'] * england['ratio']
        england['C_20'] = england['C_20'] * england['ratio']
        england['C_25'] = england['C_25'] * england['ratio']
        england['C_30'] = england['C_30'] * england['ratio']
        england['C_40'] = england['C_40'] * england['ratio']
        england['C_50'] = england['C_50'] * england['ratio']
        england['C_60'] = england['C_60'] * england['ratio']
        england['C_70'] = england['C_70'] * england['ratio']
        england['C_75'] = england['C_75'] * england['ratio']
        england['C_80'] = england['C_80'] * england['ratio']
        england['C_90'] = england['C_90'] * england['ratio']
        england['C_975'] = england['C_975'] * england['ratio']

    england = england.drop(columns=['Meta area', 'ratio']) \
                     .set_index(['area', 'day'])

    scotland = scotland_map.merge(Cpredictions.reset_index(level=1),
                                  left_on=['NHS Scotland Health Board'],
                                  right_on=['area'],
                                  how='left')

    if divide_predictions_by_population_ratio:
        scotland['C_025'] = scotland['C_025'] * scotland['ratio']
        scotland['C_10'] = scotland['C_10'] * scotland['ratio']
        scotland['C_20'] = scotland['C_20'] * scotland['ratio']
        scotland['C_25'] = scotland['C_25'] * scotland['ratio']
        scotland['C_30'] = scotland['C_30'] * scotland['ratio']
        scotland['C_40'] = scotland['C_40'] * scotland['ratio']
        scotland['C_50'] = scotland['C_50'] * scotland['ratio']
        scotland['C_60'] = scotland['C_60'] * scotland['ratio']
        scotland['C_70'] = scotland['C_70'] * scotland['ratio']
        scotland['C_75'] = scotland['C_75'] * scotland['ratio']
        scotland['C_80'] = scotland['C_80'] * scotland['ratio']
        scotland['C_90'] = scotland['C_90'] * scotland['ratio']
        scotland['C_975'] = scotland['C_975'] * scotland['ratio']

    scotland = scotland.drop(columns=['NHS Scotland Health Board', 'ratio']) \
                       .set_index(['area', 'day'])

    # We only keep the case predictions for `Highland` as a local area, not
    # Cpredictions = Cpredictions.drop(['Highland'])
    Cpredictions = Cpredictions.rename(
        index={'Highland': 'Highland (NHS Scotland Health Board)'})

    Cpredictions = Cpredictions.append(england) \
                             .append(scotland)

    # Case projections 
    england = england_map.merge(Cprojections.reset_index(level=1),
                                left_on=['Meta area'],
                                right_on=['area'],
                                how='left')

    if divide_predictions_by_population_ratio:
        england['C_025'] = england['C_025'] * england['ratio']
        england['C_10'] = england['C_10'] * england['ratio']
        england['C_20'] = england['C_20'] * england['ratio']
        england['C_25'] = england['C_25'] * england['ratio']
        england['C_30'] = england['C_30'] * england['ratio']
        england['C_40'] = england['C_40'] * england['ratio']
        england['C_50'] = england['C_50'] * england['ratio']
        england['C_60'] = england['C_60'] * england['ratio']
        england['C_70'] = england['C_70'] * england['ratio']
        england['C_75'] = england['C_75'] * england['ratio']
        england['C_80'] = england['C_80'] * england['ratio']
        england['C_90'] = england['C_90'] * england['ratio']
        england['C_975'] = england['C_975'] * england['ratio']


    england = england.drop(columns=['Meta area', 'ratio']) \
                     .set_index(['area', 'day'])

    scotland = scotland_map.merge(Cprojections.reset_index(level=1),
                                  left_on=['NHS Scotland Health Board'],
                                  right_on=['area'],
                                  how='left')

    if divide_predictions_by_population_ratio:
        scotland['C_025'] = scotland['C_025'] * scotland['ratio']
        scotland['C_10'] = scotland['C_10'] * scotland['ratio']
        scotland['C_20'] = scotland['C_20'] * scotland['ratio']
        scotland['C_25'] = scotland['C_25'] * scotland['ratio']
        scotland['C_30'] = scotland['C_30'] * scotland['ratio']
        scotland['C_40'] = scotland['C_40'] * scotland['ratio']
        scotland['C_50'] = scotland['C_50'] * scotland['ratio']
        scotland['C_60'] = scotland['C_60'] * scotland['ratio']
        scotland['C_70'] = scotland['C_70'] * scotland['ratio']
        scotland['C_75'] = scotland['C_75'] * scotland['ratio']
        scotland['C_80'] = scotland['C_80'] * scotland['ratio']
        scotland['C_90'] = scotland['C_90'] * scotland['ratio']
        scotland['C_975'] = scotland['C_975'] * scotland['ratio']


    scotland = scotland.drop(columns=['NHS Scotland Health Board', 'ratio']) \
                       .set_index(['area', 'day'])

    # We only keep the case projections for `Highland` as a local area, not
    # projections = projections.drop(['Highland'])
    Cprojections = Cprojections.rename(
        index={'Highland': 'Highland (NHS Scotland Health Board)'})

    Cprojections = Cprojections.append(england) \
                             .append(scotland)


    # Infection predictions 
    england = england_map.merge(Xpredictions.reset_index(level=1),
                                left_on=['Meta area'],
                                right_on=['area'],
                                how='left')

    if divide_predictions_by_population_ratio:
        england['C_025'] = england['C_025'] * england['ratio']
        england['C_10'] = england['C_10'] * england['ratio']
        england['C_20'] = england['C_20'] * england['ratio']
        england['C_25'] = england['C_25'] * england['ratio']
        england['C_30'] = england['C_30'] * england['ratio']
        england['C_40'] = england['C_40'] * england['ratio']
        england['C_50'] = england['C_50'] * england['ratio']
        england['C_60'] = england['C_60'] * england['ratio']
        england['C_70'] = england['C_70'] * england['ratio']
        england['C_75'] = england['C_75'] * england['ratio']
        england['C_80'] = england['C_80'] * england['ratio']
        england['C_90'] = england['C_90'] * england['ratio']
        england['C_975'] = england['C_975'] * england['ratio']

    england = england.drop(columns=['Meta area', 'ratio']) \
                     .set_index(['area', 'day'])

    scotland = scotland_map.merge(Xpredictions.reset_index(level=1),
                                  left_on=['NHS Scotland Health Board'],
                                  right_on=['area'],
                                  how='left')

    if divide_predictions_by_population_ratio:
        scotland['C_025'] = scotland['C_025'] * scotland['ratio']
        scotland['C_10'] = scotland['C_10'] * scotland['ratio']
        scotland['C_20'] = scotland['C_20'] * scotland['ratio']
        scotland['C_25'] = scotland['C_25'] * scotland['ratio']
        scotland['C_30'] = scotland['C_30'] * scotland['ratio']
        scotland['C_40'] = scotland['C_40'] * scotland['ratio']
        scotland['C_50'] = scotland['C_50'] * scotland['ratio']
        scotland['C_60'] = scotland['C_60'] * scotland['ratio']
        scotland['C_70'] = scotland['C_70'] * scotland['ratio']
        scotland['C_75'] = scotland['C_75'] * scotland['ratio']
        scotland['C_80'] = scotland['C_80'] * scotland['ratio']
        scotland['C_90'] = scotland['C_90'] * scotland['ratio']
        scotland['C_975'] = scotland['C_975'] * scotland['ratio']

    scotland = scotland.drop(columns=['NHS Scotland Health Board', 'ratio']) \
                       .set_index(['area', 'day'])

    # We only keep the case predictions for `Highland` as a local area, not
    # Xpredictions = Xpredictions.drop(['Highland'])
    Xpredictions = Xpredictions.rename(
        index={'Highland': 'Highland (NHS Scotland Health Board)'})

    Xpredictions = Xpredictions.append(england) \
                             .append(scotland)

    # Infection projections 
    england = england_map.merge(Xprojections.reset_index(level=1),
                                left_on=['Meta area'],
                                right_on=['area'],
                                how='left')

    if divide_predictions_by_population_ratio:
        england['X_025'] = england['X_025'] * england['ratio']
        england['X_10'] = england['X_10'] * england['ratio']
        england['X_20'] = england['X_20'] * england['ratio']
        england['X_25'] = england['X_25'] * england['ratio']
        england['X_30'] = england['X_30'] * england['ratio']
        england['X_40'] = england['X_40'] * england['ratio']
        england['X_50'] = england['X_50'] * england['ratio']
        england['X_60'] = england['X_60'] * england['ratio']
        england['X_70'] = england['X_70'] * england['ratio']
        england['X_75'] = england['X_75'] * england['ratio']
        england['X_80'] = england['X_80'] * england['ratio']
        england['X_90'] = england['X_90'] * england['ratio']
        england['X_975'] = england['X_975'] * england['ratio']


    england = england.drop(columns=['Meta area', 'ratio']) \
                     .set_index(['area', 'day'])

    scotland = scotland_map.merge(Xprojections.reset_index(level=1),
                                  left_on=['NHS Scotland Health Board'],
                                  right_on=['area'],
                                  how='left')

    if divide_predictions_by_population_ratio:
        scotland['X_025'] = scotland['X_025'] * scotland['ratio']
        scotland['X_10'] = scotland['X_10'] * scotland['ratio']
        scotland['X_20'] = scotland['X_20'] * scotland['ratio']
        scotland['X_25'] = scotland['X_25'] * scotland['ratio']
        scotland['X_30'] = scotland['X_30'] * scotland['ratio']
        scotland['X_40'] = scotland['X_40'] * scotland['ratio']
        scotland['X_50'] = scotland['X_50'] * scotland['ratio']
        scotland['X_60'] = scotland['X_60'] * scotland['ratio']
        scotland['X_70'] = scotland['X_70'] * scotland['ratio']
        scotland['X_75'] = scotland['X_75'] * scotland['ratio']
        scotland['X_80'] = scotland['X_80'] * scotland['ratio']
        scotland['X_90'] = scotland['X_90'] * scotland['ratio']
        scotland['X_975'] = scotland['X_975'] * scotland['ratio']


    scotland = scotland.drop(columns=['NHS Scotland Health Board', 'ratio']) \
                       .set_index(['area', 'day'])

    # We only keep the case projections for `Highland` as a local area, not
    # projections = projections.drop(['Highland'])
    if divide_predictions_by_population_ratio:
        Xprojections = Xprojections.rename(
            index={'Highland': 'Highland (NHS Scotland Health Board)'})

    Xprojections = Xprojections.append(england) \
                             .append(scotland)

    # The original reproduction numbers and case predictions data frames might
    # already have been reinflated, so that we append duplicate rows. Duplicate
    # rows are removed here. Note that we remove duplicates by index, and not
    # by the entire row, as that requires floating point precision equality
    # (which we don't get from reading from CSV).
    reproduction_numbers = reproduction_numbers.loc[
        ~reproduction_numbers.index.duplicated(keep='first')]
    exceedance_probs = exceedance_probs.loc[
        ~exceedance_probs.index.duplicated(keep='first')]
    weekly = weekly.loc[
        ~weekly.index.duplicated(keep='first')]
    Cpredictions = Cpredictions.loc[
        ~Cpredictions.index.duplicated(keep='first')]
    Cprojections = Cprojections.loc[
        ~Cprojections.index.duplicated(keep='first')]
    Xpredictions = Xpredictions.loc[
        ~Xpredictions.index.duplicated(keep='first')]
    Xprojections = Xprojections.loc[
        ~Xprojections.index.duplicated(keep='first')]

    return reproduction_numbers.sort_index().dropna(), \
           exceedance_probs.sort_index().dropna(), \
           weekly.sort_index().dropna(), \
           Cpredictions.sort_index().dropna(), \
           Cprojections.sort_index().dropna(), \
           Xpredictions.sort_index().dropna(), \
           Xprojections.sort_index().dropna() 


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) == 14:
        process_csvs(input_r_filename=args[0],
                     input_exceed_filename=args[1],
                     input_weekly_filename=args[2],
                     input_Cpred_filename=args[3],
                     input_Cproj_filename=args[4],
                     input_Xpred_filename=args[5],
                     input_Xproj_filename=args[6],
                     output_r_filename=args[7],
                     output_exceed_filename=args[8],
                     output_weekly_filename=args[9],
                     output_Cpred_filename=args[10],
                     output_Cproj_filename=args[11],
                     output_Xpred_filename=args[12],
                     output_Xproj_filename=args[13])
    elif len(args) == 7:
        process_csvs(input_r_filename=args[0],
                     input_exceed_filename=args[1],
                     input_weekly_filename=args[2],
                     input_Cpred_filename=args[3],
                     input_Cproj_filename=args[4],
                     input_Xpred_filename=args[5],
                     input_Xproj_filename=args[6])
    else:
        raise Error("Number of arguments has to be 7 or 14")
