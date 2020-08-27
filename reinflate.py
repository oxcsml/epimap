import pandas as pd
import sys


def process_csvs(input_r_filename: str='fits/0_Rt.csv',
                 input_pred_filename: str='fits/0_Cpred.csv',
                 input_proj_filename: str='fits/0_Cproj.csv',
                 output_r_filename: str='website/Rt.csv',
                 output_pred_filename: str='website/Cpred.csv',
                 output_proj_filename: str='website/Cproj.csv'
):

    reproduction_numbers = pd.read_csv(input_r_filename) \
                             .rename(columns={'Date': 'day'})
    reproduction_numbers['day'] = pd.to_datetime(reproduction_numbers['day'],
                                                 format='%Y-%m-%d')
    reproduction_numbers = reproduction_numbers.set_index(['area', 'day'])

    predictions = pd.read_csv(input_pred_filename) \
                         .rename(columns={'Date': 'day'})
    predictions['day'] = pd.to_datetime(predictions['day'],
                                             format='%Y-%m-%d')
    predictions = predictions.set_index(['area', 'day'])


    projections = pd.read_csv(input_proj_filename) \
                         .rename(columns={'Date': 'day'})
    projections['day'] = pd.to_datetime(projections['day'],
                                             format='%Y-%m-%d')
    projections = projections.set_index(['area', 'day'])

    reproduction_numbers, predictions, projections = reinflate(
        reproduction_numbers, predictions, projections)

    reproduction_numbers.index.names = ['area', 'Date']
    predictions.index.names = ['area', 'Date']
    projections.index.names = ['area', 'Date']

    reproduction_numbers.to_csv(output_r_filename, float_format='%.5f')
    predictions.to_csv(output_pred_filename, float_format='%.5f')
    projections.to_csv(output_proj_filename, float_format='%.5f')


def reinflate(reproduction_numbers: pd.DataFrame,
              predictions: pd.DataFrame,
              projections: pd.DataFrame,
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
    england = england_map.merge(predictions.reset_index(level=1),
                                left_on=['Meta area'],
                                right_on=['area'],
                                how='left')

    if divide_predictions_by_population_ratio:
        england['C_025'] = england['C_025'] * england['ratio']
        england['C_25']  = england['C_25']  * england['ratio']
        england['C_50']  = england['C_50']  * england['ratio']
        england['C_75']  = england['C_75']  * england['ratio']
        england['C_975'] = england['C_975'] * england['ratio']

    england = england.drop(columns=['Meta area', 'ratio']) \
                     .set_index(['area', 'day'])

    scotland = scotland_map.merge(predictions.reset_index(level=1),
                                  left_on=['NHS Scotland Health Board'],
                                  right_on=['area'],
                                  how='left')

    if divide_predictions_by_population_ratio:
        scotland['C_025'] = scotland['C_025'] * scotland['ratio']
        scotland['C_25']  = scotland['C_25']  * scotland['ratio']
        scotland['C_50']  = scotland['C_50']  * scotland['ratio']
        scotland['C_75']  = scotland['C_75']  * scotland['ratio']
        scotland['C_975'] = scotland['C_975'] * scotland['ratio']

    scotland = scotland.drop(columns=['NHS Scotland Health Board', 'ratio']) \
                       .set_index(['area', 'day'])

    # We only keep the case predictions for `Highland` as a local area, not
    # predictions = predictions.drop(['Highland'])
    if divide_predictions_by_population_ratio:
        predictions = predictions.rename(
            index={'Highland': 'Highland (NHS Scotland Health Board)'})

    predictions = predictions.append(england) \
                             .append(scotland)

    # Case projections 
    england = england_map.merge(projections.reset_index(level=1),
                                left_on=['Meta area'],
                                right_on=['area'],
                                how='left')

    if divide_predictions_by_population_ratio:
        england['C_025'] = england['C_025'] * england['ratio']
        england['C_25']  = england['C_25']  * england['ratio']
        england['C_50']  = england['C_50']  * england['ratio']
        england['C_75']  = england['C_75']  * england['ratio']
        england['C_975'] = england['C_975'] * england['ratio']

    england = england.drop(columns=['Meta area', 'ratio']) \
                     .set_index(['area', 'day'])

    scotland = scotland_map.merge(projections.reset_index(level=1),
                                  left_on=['NHS Scotland Health Board'],
                                  right_on=['area'],
                                  how='left')

    if divide_predictions_by_population_ratio:
        scotland['C_025'] = scotland['C_025'] * scotland['ratio']
        scotland['C_25']  = scotland['C_25']  * scotland['ratio']
        scotland['C_50']  = scotland['C_50']  * scotland['ratio']
        scotland['C_75']  = scotland['C_75']  * scotland['ratio']
        scotland['C_975'] = scotland['C_975'] * scotland['ratio']

    scotland = scotland.drop(columns=['NHS Scotland Health Board', 'ratio']) \
                       .set_index(['area', 'day'])

    # We only keep the case projections for `Highland` as a local area, not
    # projections = projections.drop(['Highland'])
    if divide_predictions_by_population_ratio:
        projections = projections.rename(
            index={'Highland': 'Highland (NHS Scotland Health Board)'})

    projections = projections.append(england) \
                             .append(scotland)


    # The original reproduction numbers and case predictions data frames might
    # already have been reinflated, so that we append duplicate rows. Duplicate
    # rows are removed here. Note that we remove duplicates by index, and not
    # by the entire row, as that requires floating point precision equality
    # (which we don't get from reading from CSV).
    reproduction_numbers = reproduction_numbers.loc[
        ~reproduction_numbers.index.duplicated(keep='first')]
    predictions = predictions.loc[
        ~predictions.index.duplicated(keep='first')]
    projections = projections.loc[
        ~projections.index.duplicated(keep='first')]

    return reproduction_numbers.sort_index(), \
           predictions.sort_index(), \
           projections.sort_index()


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) == 6:
        process_csvs(input_r_filename=args[0],
                     input_pred_filename=args[1],
                     input_proj_filename=args[2],
                     output_r_filename=args[3],
                     output_pred_filename=args[4],
                     output_proj_filename=args[5]
        )
    elif len(args) == 3:
        process_csvs(input_r_filename=args[0],
                     input_pred_filename=args[1],
                     input_proj_filename=args[2]
        )
    else:
        process_csvs()
