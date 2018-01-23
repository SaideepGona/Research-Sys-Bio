#!/usr/bin/env python3
from argparse import ArgumentParser
from datetime import date
from fnmatch import fnmatch
from ftplib import FTP
from io import BytesIO
import json
from pathlib import Path, PurePosixPath
import pickle
from urllib.parse import urlunparse, urlparse
import tarfile
from typing import Any, Dict, List, Optional, Tuple

from data_path_utils import create_data_path
from lxml import etree as ET

from query_download_utils import create_ftp_session, query_all_accession_numbers
from utils import DOWNLOAD_PATH, first, strip_prefix

METADATA_PATH = DOWNLOAD_PATH / 'metadata'

SRA_PATTERN = '*.sra'

NAMESPACE_KEY = 'miniml'

SAMPLE_TAG = f'{NAMESPACE_KEY}:Sample'
SERIES_TAG = f'{NAMESPACE_KEY}:Series'
CHANNEL_TAG = f'{NAMESPACE_KEY}:Channel'
SOURCE_TAG = f'{NAMESPACE_KEY}:Source'
ORGANISM_TAG = f'{NAMESPACE_KEY}:Organism'
CHARACTERISTICS_TAG = f'{NAMESPACE_KEY}:Characteristics'
SUPPLEMENTARY_DATA_TAG = f'{NAMESPACE_KEY}:Supplementary-Data'
PUBMED_ID_TAG = f'{NAMESPACE_KEY}:Pubmed-ID'
LIBRARY_STRATEGY_TAG = f'{NAMESPACE_KEY}:Library-Strategy'

SOURCE_XPATH = f'{CHANNEL_TAG}/{SOURCE_TAG}[1]/text()'
ORGANISM_XPATH = f'{CHANNEL_TAG}/{ORGANISM_TAG}[1]/text()'
CELL_TYPE_XPATH = f'{CHANNEL_TAG}/{CHARACTERISTICS_TAG}[@tag="cell type"]/text()'
CELL_ID_XPATH = f'{CHANNEL_TAG}/{CHARACTERISTICS_TAG}[@tag="cell_id"]/text()'
SRA_FTP_XPATH = f'{SUPPLEMENTARY_DATA_TAG}[@type="SRA Experiment"]/text()'
PUBMED_ID_XPATH = f'{SERIES_TAG}/{PUBMED_ID_TAG}/text()'
LIBRARY_STRATEGY_XPATH = f'{LIBRARY_STRATEGY_TAG}/text()'

# TODO: Figure out a nicer way to build XPath queries from tag names; this is a huge mess

SUBMISSION_DATE_XPATH = f'{NAMESPACE_KEY}:Status/{NAMESPACE_KEY}:Submission-Date[1]/text()'
RELEASE_DATE_XPATH = f'{NAMESPACE_KEY}:Status/{NAMESPACE_KEY}:Release-Date[1]/text()'
LAST_UPDATE_DATE_XPATH = f'{NAMESPACE_KEY}:Status/{NAMESPACE_KEY}:Last-Update-Date[1]/text()'

REQUIRED_ORGANISM = 'Mus musculus'

class SampleInfo:
    series: str
    iid: str
    source: str
    cell_type: Optional[str]
    cell_id: Optional[str]
    submission_date: date
    release_date: date
    last_update_date: date
    sra_path: str
    pubmed_id: str

    def __init__(self):
        for name in self.__annotations__:
            setattr(self, name, None)

    def to_dict(self) -> Dict[str, Any]:
        return {name: str(getattr(self, name)) for name in self.__annotations__}

def get_local_family_file_path(accn: str) -> Path:
    local_family_file_path = METADATA_PATH / accn / f'{accn}_family.xml'
    local_family_file_path.parent.mkdir(parents=True, exist_ok=True)
    return local_family_file_path

def download_family_file(accn: str, ftp_session: Optional[FTP]=None):
    accn_prefix = int(strip_prefix(accn, 'GSE')) // 1000
    miniml_path = f'/geo/series/GSE{accn_prefix}nnn/{accn}/miniml'
    if ftp_session:
        ftp = ftp_session
    else:
        print('Connecting to FTP server')
        ftp = create_ftp_session()
    print(f'Changing FTP directory to {miniml_path}')
    ftp.cwd(miniml_path)
    b = BytesIO()
    ftp.retrbinary(f'RETR {accn}_family.xml.tgz', b.write)
    b.seek(0)

    if ftp_session is None:
        # We created the session inside this function call, so log out
        ftp.quit()

    local_family_file_path = get_local_family_file_path(accn)

    # unzip tgz file to family.xml
    tar = tarfile.open(fileobj=b, mode='r:gz')
    for member in tar.getmembers():
        # only grab the family.xml file from the bundle
        if member.name == local_family_file_path.name:
            break
    else:
        raise ValueError(f'No {local_family_file_path.name} file in archive')

    f = tar.extractfile(member)
    family_file = f.read()

    with open(local_family_file_path, 'wb') as o:
        o.write(family_file)

def get_first_value(node, xpath: str, nsmap: dict) -> Optional[str]:
    values = node.xpath(xpath, namespaces=nsmap)
    if values:
        return first(values).strip()

def parse_family_file_contents(family_file_contents: bytes, accn: str) -> List[SampleInfo]:
    sample_data: List[SampleInfo] = []

    try:
        famRoot = ET.fromstring(family_file_contents)
    except Exception:
        print('XML parse error:')
        print(family_file_contents[:200])
        return []
    nsmap = {NAMESPACE_KEY: famRoot.nsmap[None]}

    pubmed_ids = famRoot.xpath(PUBMED_ID_XPATH, namespaces=nsmap)
    if pubmed_ids:
        pubmed_id = pubmed_ids[0]
    else:
        pubmed_id = None

    for sample in famRoot.xpath(SAMPLE_TAG, namespaces=nsmap):
        s = SampleInfo()
        s.series = accn
        s.iid = sample.attrib['iid']

        try:
            # Need both of these to be present; `sample.xpath(...)` will return an empty list
            # if not, and our `[0]` access will throw an IndexError
            s.source = sample.xpath(SOURCE_XPATH, namespaces=nsmap)[0]
            organism = sample.xpath(ORGANISM_XPATH, namespaces=nsmap)[0]
        except IndexError:
            print('No Source or organism found')
            continue

        if organism != REQUIRED_ORGANISM:
            print(f'invalid organism ({organism}) found')
            continue

        s.cell_type = get_first_value(sample, CELL_TYPE_XPATH, nsmap)
        s.cell_id = get_first_value(sample, CELL_ID_XPATH, nsmap)
        s.submission_date = get_first_value(sample, SUBMISSION_DATE_XPATH, nsmap)
        s.release_date = get_first_value(sample, RELEASE_DATE_XPATH, nsmap)
        s.last_update_date = get_first_value(sample, LAST_UPDATE_DATE_XPATH, nsmap)
        s.sra_path = get_first_value(sample, SRA_FTP_XPATH, nsmap)
        s.pubmed_id = pubmed_id

        sample_data.append(s)

    return sample_data

def find_sra_file(si: SampleInfo, ftp: FTP) -> Optional[str]:
    print('Finding SRA file for', si.iid)
    if si.sra_path is None:
        return None

    ftp_url_pieces = urlparse(si.sra_path)
    try:
        # It's a slight pain to wrap these Path objects in `str` for use
        # in the FTP module, but the path manipulation still makes it
        # worth it in my opinion
        path = PurePosixPath(ftp_url_pieces.path)
        ftp.cwd(str(path))
        dirs = list(ftp.mlsd())
        # '.', '..', directory name we care about
        assert len(dirs) == 3
        name, attrs = dirs[-1]
        # '.' and '..' should always be first, but let's make sure just in case
        assert name.strip('.')
        assert attrs['type'] == 'dir'
        new_path = path / name
        ftp.cwd(str(new_path))
        files = list(ftp.mlsd())
        assert len(files) == 3
        filename, file_attrs = files[-1]
        assert file_attrs['type'] == 'file'
        assert fnmatch(filename, SRA_PATTERN)

        url_piece_list = list(ftp_url_pieces)
        url_piece_list[2] = str(new_path / filename)

        return urlunparse(url_piece_list)
    except Exception:
        return None

def parse_family_file(accn: str, ftp: Optional[FTP]=None, sra_paths: bool=False) -> List[SampleInfo]:
    """
    :param accn: Series ID, like 'GSE79578'
    :return: List of SampleInfo objects
    """
    local_family_file_path = get_local_family_file_path(accn)
    if not local_family_file_path.is_file():
        download_family_file(accn, ftp)

    with open(local_family_file_path, 'rb') as f:
        family_file_contents = f.read()

    contents = parse_family_file_contents(family_file_contents, accn)
    if sra_paths:
        # SRA paths are directories here; get actual files
        for sample_info in contents:
            sample_info.sra_path = find_sra_file(sample_info, ftp)

    return contents

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument(
        'accession_number',
        nargs='*',
        help='Example: GSE54013. If not specified, all available accession numbers will be queried.',
    )
    p.add_argument('--sra-paths', action='store_true')
    args = p.parse_args()

    accession_numbers = args.accession_number or query_all_accession_numbers(limit=None)

    data_path = create_data_path('parse_family_file')

    with create_ftp_session() as ftp:
        accn_data: Dict[str, List[SampleInfo]] = {}
        for accn in accession_numbers:
            accn_data[accn] = parse_family_file(accn, ftp, args.sra_paths)

    accn_data_adj = {
        accn: [sample_info.to_dict() for sample_info in data]
        for accn, data in accn_data.items()
    }

    json_data_path = data_path / 'metadata.json'
    print('Saving metadata to', json_data_path)
    with open(json_data_path, 'w') as f:
        json.dump(accn_data_adj, f)

    pickle_data_path = data_path / 'metadata.pickle'
    print('Saving full metadata to', pickle_data_path)
    with open(pickle_data_path, 'wb') as f:
        pickle.dump(accn_data, f)
