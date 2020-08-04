import pandas as pd
import sys


def process_csvs(input_r_filename: str='fits/0_Rt.csv',
                 input_cases_filename: str='fits/0_Cproj.csv',
                 output_r_filename: str='website/Rt.csv',
                 output_cases_filename: str='website/Cproj.csv'):

    reproduction_numbers = pd.read_csv(input_r_filename) \
                             .set_index('area')
    case_predictions = pd.read_csv(input_cases_filename) \
                         .rename(columns={'Date': 'day'})
    case_predictions['day'] = pd.to_datetime(case_predictions['day'],
                                             format='%Y-%m-%d')
    case_predictions = case_predictions.set_index(['area', 'day'])

    reproduction_numbers, case_predictions = reinflate(
        reproduction_numbers, case_predictions)

    case_predictions.index.names = ['area', 'Date']

    reproduction_numbers.to_csv(output_r_filename, float_format='%.5f')
    case_predictions.to_csv(output_cases_filename, float_format='%.5f')


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
