#!/usr/bin/env python3
from argparse import ArgumentParser

from data_path_utils import create_data_path
import pandas as pd

from parse_family_file import parse_family_file
from query_download_utils import create_ftp_session, query_all_accession_numbers

def main(sra_paths: bool):
    data_path = create_data_path('bulk_family_file_query')

    series_data = []
    with create_ftp_session() as ftp:
        for accn in query_all_accession_numbers():
            series_data.extend(parse_family_file(accn, ftp=ftp, sra_paths=sra_paths))

    d = pd.DataFrame([si.to_dict() for si in series_data])
    series_data_path = data_path / 'series_data.csv'
    print('Saving series data to', series_data_path)
    d.to_csv(series_data_path)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('--sra-paths', action='store_true')
    args = p.parse_args()

    main(args.sra_paths)
