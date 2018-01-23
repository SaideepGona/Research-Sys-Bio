from ftplib import FTP
from functools import lru_cache
from typing import Iterable, List, Optional, Tuple, Union
import urllib.error
import urllib.parse
import urllib.request

from lxml import etree as ET

FTP_HOST = 'ftp.ncbi.nlm.nih.gov'
FTP_PATH_TEMPLATE = '/geo/series/GSE{accn_prefix}nnn/GSE{accn}/suppl'

def create_ftp_session(host=FTP_HOST) -> FTP:
    ftp = FTP(host)
    ftp.set_pasv(True)
    # anonymous
    ftp.login()
    return ftp

REQUEST_TEMPLATE_PIECES = (
    'https',
    'eutils.ncbi.nlm.nih.gov',
    '/entrez/eutils/{action}.fcgi',
    '',
    '{query}',
    '',
)
ACTION_INDEX = 2
QUERY_INDEX = 4

UID_START = 200000000

# The default limit is 20 if a limit isn't specified in the query -- so we need
# to override with something large. Name of the variable means "the value we use
# for an unlimited query". Queries are limited server-side to 100,000 as per
# https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.retmax and we will likely
# not ever encounter more expression datasets than this, so let's just use that.
REQUEST_LIMIT = 100_000

def build_request(query_str: str, action: str) -> str:
    """
    Constructs a HTTP(S) request URL from `REQUEST_TEMPLATE_PIECES`, replacing
    the `query` piece with the provided `query_str`.
    """
    request_pieces = list(REQUEST_TEMPLATE_PIECES)
    request_pieces[ACTION_INDEX] = request_pieces[ACTION_INDEX].format(action=action)
    request_pieces[QUERY_INDEX] = query_str
    return urllib.parse.urlunparse(request_pieces)

class QueryNode:
    op: str
    terms: List[Union['QueryNode', str]]

    def __init__(self, *terms):
        # Want to be able to alter this after creation; argument unpacking
        # provides a tuple
        self.terms = list(terms)

    def build_query_str(self):
        # Add some spaces around operator, to use in joining terms
        sep = f' {self.op} '
        term_strs = [
            term if isinstance(term, str) else term.build_query_str()
            for term in self.terms
        ]
        query_str = '({})'.format(sep.join(term_strs))
        return query_str

class AndNode(QueryNode):
    op = 'AND'

class OrNode(QueryNode):
    op = 'OR'

def build_query_str():
    q = AndNode(
        'mouse[orgn]',
        OrNode(
            'single cell[all]',
            'single-cell[all]',
            'scRNA-seq[all]',
        ),
        'RNA[all]',
        'expression profiling by high throughput sequencing[GTYP]',
    )
    return q.build_query_str()

def query_all_accession_numbers(limit: Optional[int]=REQUEST_LIMIT) -> Iterable[str]:
    # build query
    limit_value = limit or REQUEST_LIMIT
    query_pieces = [
        ('db', 'gds'),
        ('term', build_query_str()),
        ('retmax', limit_value),
    ]
    query_str = urllib.parse.urlencode(query_pieces)
    search = build_request(query_str, 'esearch')

    print('Querying', search)
    f = urllib.request.urlopen(search)
    contents = f.read()
    root = ET.fromstring(contents)

    for id_node in root.find('IdList').findall('Id'):
        try:
            uid = int(id_node.text)
        except ValueError:
            print(f"Couldn't parse UID {id_node.text}")
            continue
        # All records have UID 2000xxxxx, the last 5 digits of which correspond to the accession number
        if uid < UID_START:
            print(f'Unrecognized UID: {uid}')
            continue

        accn = uid - UID_START

        yield f'GSE{accn}'

@lru_cache(maxsize=None)
def query_srr_id(srr_id: str) -> str:
    """
    :param srr_id: filename base of a read archive, e.g. "SRR3319480"
    :return: Internal GEO ID for that sample, which is a numeric string,
    like "2399278" for that SRR number

    Intended only for use by query_srr_series_sample.
    """
    query_pieces = [
        ('db', 'sra'),
        ('term', srr_id),
    ]
    query_str = urllib.parse.urlencode(query_pieces)
    search = build_request(query_str, 'esearch')

    print('Querying', search)
    f = urllib.request.urlopen(search)
    contents = f.read()
    root = ET.fromstring(contents)

    id_nodes = root.xpath('//Id')
    assert len(id_nodes) == 1
    return id_nodes[0].text

@lru_cache(maxsize=None)
def query_srr_series_sample(srr_id: str) -> Tuple[str, str]:
    """
    :param srr_id: filename base of a read archive, e.g. "SRR3319480"
    :return: 2-tuple:
     [0] Series number, e.g. "GSE79812"
     [1] Sample number, e.g. "GSM2104026"
    """
    geo_internal_id = query_srr_id(srr_id)
    query_pieces = [
        ('db', 'sra'),
        ('id', geo_internal_id),
    ]
    query_str = urllib.parse.urlencode(query_pieces)
    search = build_request(query_str, 'efetch')

    print('Querying', search)
    f = urllib.request.urlopen(search)
    contents = f.read()
    root = ET.fromstring(contents)

    study_nodes = root.xpath('//STUDY')
    assert len(study_nodes) == 1
    study_id = study_nodes[0].attrib['alias']

    sample_nodes = root.xpath('//SAMPLE')
    assert len(sample_nodes) == 1
    sample_id = sample_nodes[0].attrib['alias']

    return study_id, sample_id
