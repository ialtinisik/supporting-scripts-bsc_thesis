import pandas as pd
from pathlib import Path
import numpy as np


def format_wellID(id: str):
    ''' maps well without 0 index to with 0 '''
    if not isinstance(id, str):
        return id
    elif len(id) == 3:
        return id
    elif len(id) == 2:
        return str(id[0] + '0' + id[1])
    else:
        raise NotImplementedError('unexpected well input')


def make_wellplate_indices():
    row = "ABCDEFGH"
    col = range(1, 13)
    return [format_wellID(str(r)+str(c)) for r in row for c in col]


def offset_reactor_block_pos(plate_id, peptides_per_plate=15):
    wellplate = make_wellplate_indices()
    offset = peptides_per_plate * plate_id
    return wellplate[offset: offset + peptides_per_plate]

    # real_pos = [
    #     offset_reactor_block_pos(plate_id=i)
    #     for i in range(len(plate_file_paths))
    # ]
    # real_reactor_block_loc = [
    #     pos for well_pos in real_pos for pos in well_pos
    # ]


def map_pos_to_offset_pos(well_id, plate_id):
    if pd.isna(well_id):
        return 'Z99'
    wellplate = make_wellplate_indices()
    offset_wellplate = offset_reactor_block_pos(plate_id=plate_id)
    index = wellplate.index(well_id)
    return offset_wellplate[index]


def map_well_to_seq(seq_file_path: str, plate_file_dir: str):
    """ maps (resolves) plate reader output files to sequence to obtain
        (`sequence`, `concentration`, `CSC` (= `normintegral`) output.csv

        uses Mapping_PSM0009.csv for mapping, which are not expected to change:
            - has columns: `platereader_loc`, `reactor_block_loc`, `platereader_conc`

        params:
            seq_file_path: path to file with (sequence, reactor_block_pos) info
                - expects columns `reactor_block_loc`, `sequence`
            plate_file_dir: path to directory with plate reader outputs (.txt)
                - expects multiple .txt files
    """
    # seq_file_path = './data/PSM0009_synlocs_to_seqs.csv'
    seq_df = pd.read_csv(seq_file_path)
    print('seq_df', seq_df.columns)
    seq_df['reactor_block_loc'] = seq_df['reactor_block_loc'].apply(
        format_wellID
    )

    # map_df and map_file_path not expected to change across experiments
    map_df = pd.read_csv('./data/reactor_to_well_mapping.csv')
    map_df['reactor_block_loc'] = map_df['reactor_block_loc'].apply(
        format_wellID)
    map_df['platereader_loc'] = map_df['platereader_loc'].apply(
        format_wellID)

    # print('map_df', map_df.columns, '\n', 'map_df', map_df.head())

    # TODO: add analysis script here to get ['Well', 'CSC' or 'normintegral']
    # or easier: read analysed file from CSCanyalsis script directly and
    # only keep `platereader_loc` and `csc/integral/normintegral` columns

    # plate_file_dir: directory with plate reader .txt nephelometry files
    # plate_file_dir = './data/Measurements_PSM0009'

    # TODO verify mergeing does not break colnames through inner join
    # check assumption that these are in order

    plate_file_paths = list(Path(plate_file_dir).glob('*.txt'))
    plate_df = pd.concat(
        [pd.read_csv(path, sep='\t', header=[0, 1])
         for path in plate_file_paths],
        axis=0).reset_index(drop=True)

    plate_df['plate_id'] = [i for i in range(len(plate_file_paths))
                            for _ in range(96)]

    plate_df['plate_id'] = [i for i in range(len(plate_file_paths))
                            for _ in range(96)]

    # plate_df['real_reactor_block_loc'] = plate_df['reactor_block_loc'].apply(
    #     map_pos_to_offset_pos)

    # print(real_reactor_block_loc)

    # plate_df = plate_df.dropna([)

    # flatten the multi-index & rename (from plates.txt with 2 header rows)
    plate_df.columns = [' '.join(col).strip().replace('Raw Data ', '')
                        for col in plate_df.columns.values]
    plate_df = plate_df.rename(
        {'Well Unnamed: 0_level_1': 'platereader_loc'}, axis=1)

    plate_df.to_csv('./data/debug_platedf.csv')

    # merge (join) measurements (plate_df) with map_df to relate
    # platereader loc to reactor_block loc
    # TODO check merge INNER is fine due to 4x96 measurements
    joined_df = pd.merge(plate_df, map_df, on='platereader_loc', how='left')

    # remap reactor block pos based on plate_id
    joined_df['reactor_block_loc'] = joined_df.apply(
        lambda row: map_pos_to_offset_pos(
            row['reactor_block_loc'], row['plate_id']), axis=1)

    # merge (join) this with seq_df to relate reactor_block_loc to sequence
    joined_df = pd.merge(
        joined_df, seq_df, on='reactor_block_loc', how='left')

    # print('plate_df columns:', plate_df.columns.values, '\n',
    #       'seq_df columns:', seq_df.columns.values, '\n',
    #       'map_df columns:', map_df.columns.values, '\n',
    #       'merged_df columns:', joined_df.columns.values)

    # optional: remove excess columns
    # joined_df = joined_df.loc[:, ['sequence',
    #                               'platereader_conc', 'reactor_block_loc']]
    # , 'normintegral', 'integral']]
    print(joined_df.sample(5))

    return joined_df


if __name__ == "__main__":
    experiment_name = 'PSM0009'
    merged_df = map_well_to_seq(
        seq_file_path=f'./data/{experiment_name}_synlocs_to_seqs.csv',
        plate_file_dir=f'./data/Measurements_{experiment_name}',
    )

    merged_df.to_csv(f'./data/seq_exp_{experiment_name}.csv')
