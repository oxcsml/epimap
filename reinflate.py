import pandas as pd
import sys


def process_csvs(input_r_filename: str='fits/0_Rt.csv',
                 input_cases_filename: str='fits/0_Cproj.csv',
                 output_r_filename: str='website/Rt.csv',
                 output_cases_filename: str='website/Cproj.csv'):

    reproduction_numbers = pd.read_csv(input_r_filename) \
                             .rename(columns={'Date': 'day'})
    reproduction_numbers['day'] = pd.to_datetime(reproduction_numbers['day'],
                                                 format='%Y-%m-%d')
    reproduction_numbers = reproduction_numbers.set_index(['area', 'day'])

    case_predictions = pd.read_csv(input_cases_filename) \
                         .rename(columns={'Date': 'day'})
    case_predictions['day'] = pd.to_datetime(case_predictions['day'],
                                             format='%Y-%m-%d')
    case_predictions = case_predictions.set_index(['area', 'day'])

    reproduction_numbers, case_predictions = reinflate(
        reproduction_numbers, case_predictions)

    case_predictions.index.names = ['area', 'Date']
    reproduction_numbers.index.names = ['area', 'Date']

    reproduction_numbers.to_csv(output_r_filename, float_format='%.5f')
    case_predictions.to_csv(output_cases_filename, float_format='%.5f')


def reinflate(reproduction_numbers: pd.DataFrame,
              case_predictions: pd.DataFrame,
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
    england = england_map.merge(reproduction_numbers.reset_index(level=1),
                                left_on=['Meta area'],
                                right_on=['area'],
                                how='left')

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

    # `Highland` is both the name of an NHS Scotland Health Board and an area
    # in the health board. To avoid a duplicate area key, we drop it as a
    # health board, and keep it as a local area (in this case the two rows
    # would have been equal in any event).
    if divide_predictions_by_population_ratio:
        reproduction_numbers = reproduction_numbers.rename(
            index={'Highland': 'Highland (NHS Scotland Health Board)'})

    reproduction_numbers = reproduction_numbers.append(england) \
                                               .append(scotland)

    # Case predictions
    england = england_map.merge(case_predictions.reset_index(level=1),
                                left_on=['Meta area'],
                                right_on=['area'],
                                how='left')

    if divide_predictions_by_population_ratio:
        england['C_lower'] = england['C_lower'] * england['ratio']
        england['C_median'] = england['C_median'] * england['ratio']
        england['C_upper'] = england['C_upper'] * england['ratio']

    england = england.drop(columns=['Meta area', 'ratio']) \
                     .set_index(['area', 'day'])

    scotland = scotland_map.merge(case_predictions.reset_index(level=1),
                                  left_on=['NHS Scotland Health Board'],
                                  right_on=['area'],
                                  how='left')

    if divide_predictions_by_population_ratio:
        scotland['C_lower'] = scotland['C_lower'] * scotland['ratio']
        scotland['C_median'] = scotland['C_median'] * scotland['ratio']
        scotland['C_upper'] = scotland['C_upper'] * scotland['ratio']

    scotland = scotland.drop(columns=['NHS Scotland Health Board', 'ratio']) \
                       .set_index(['area', 'day'])

    # We only keep the case predictions for `Highland` as a local area, not
    # case_predictions = case_predictions.drop(['Highland'])
    if divide_predictions_by_population_ratio:
        case_predictions = case_predictions.rename(
            index={'Highland': 'Highland (NHS Scotland Health Board)'})

    case_predictions = case_predictions.append(england) \
                                       .append(scotland)

    # The original reproduction numbers and case predictions data frames might
    # already have been reinflated, so that we append duplicate rows. Duplicate
    # rows are removed here. Note that we remove duplicates by index, and not
    # by the entire row, as that requires floating point precision equality
    # (which we don't get from reading from CSV).
    reproduction_numbers = reproduction_numbers.loc[
        ~reproduction_numbers.index.duplicated(keep='first')]
    case_predictions = case_predictions.loc[
        ~case_predictions.index.duplicated(keep='first')]

    return reproduction_numbers.sort_index(), case_predictions.sort_index()


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) == 4:
        process_csvs(input_r_filename=args[0],
                     input_cases_filename=args[1],
                     output_r_filename=args[2],
                     output_cases_filename=args[3])
    elif len(args) == 2:
        process_csvs(input_r_filename=args[0],
                     input_cases_filename=args[1])
    else:
        process_csvs()
