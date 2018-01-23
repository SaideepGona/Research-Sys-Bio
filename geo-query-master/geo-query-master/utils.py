from pathlib import Path
from typing import Iterable, TypeVar

TIMESTAMP_FORMAT = '%Y%m%d-%H%M%S'

DATA_PATH = Path('data')
OUTPUT_PATH = Path('output')
SLURM_PATH = Path('slurm')
DOWNLOAD_PATH = Path('download')
DOWNLOAD_PATH_NEW = DOWNLOAD_PATH / 'new'
DOWNLOAD_PATH_PROCESSED = DOWNLOAD_PATH / 'processed'

T = TypeVar('T')

def first(iterable: Iterable[T]) -> T:
    return next(iter(iterable))

del T

def sorted_set_op(items, func):
    sets = [set(item) for item in items]
    data = func(*sets)
    return sorted(data)

def sorted_intersection(*items):
    return sorted_set_op(items, set.intersection)

def sorted_union(*items):
    return sorted_set_op(items, set.union)

def strip_prefix(string: str, prefix: str) -> str:
    """
    :param string: String to strip `prefix` from, if present
    :param prefix: Prefix to remove from `string`
    :return:
    """
    if string.startswith(prefix):
        return string[len(prefix):]
    return string

def strip_suffix(string: str, suffix: str) -> str:
    """
    :param string: String to strip `suffix` from, if present
    :param suffix: Suffix to remove from `string`
    :return:
    """
    if string.endswith(suffix):
        return string[:-len(suffix)]
    return string
