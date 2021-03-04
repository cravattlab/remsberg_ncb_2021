from rifl.exceptions import FilterNotFoundException
import rifl.utils as utils

from collections import namedtuple
from typing import List, Literal, Callable, Union, Optional
from dataclasses import dataclass
import functools

FilterFunction = Callable[..., Union[bool, List[bool]]]
FilterLevel = Literal['individual', 'experiment', 'aggregate']


@dataclass
class Requirement:
    database_fields: Union[str, List[str]]
    dataframe_columns: Union[str, List[str]]

@dataclass
class Filter:
    fn: FilterFunction
    name: str
    requires: Requirement
    level: FilterLevel

    @classmethod
    def from_dict(cls, user_defined_filter: Union[str, dict], filters_module):
        f_kwargs = {}
        f_args = []

        if isinstance(user_defined_filter, dict):
            f, f_kwargs = next(iter(user_defined_filter.items()))

            if not isinstance(f_kwargs, dict):
                f_kwargs = {}
                f_args = f_kwargs
        else:
            f = user_defined_filter

        filter_obj = getattr(filters_module, f)

        if filter_obj is None:
            raise FilterNotFoundException
        
        filter_obj.fn = functools.partial(filter_obj.fn, *f_args, **f_kwargs)
        return filter_obj

    @property
    def dataframe_reqs(self):
        reqs = self.requires.dataframe_columns
        return reqs if isinstance(reqs, list) else [reqs]

def filter_fun(requires: Union[Requirement, str], level: Union[FilterLevel, List[FilterLevel]]):
    if isinstance(requires, str):
        requires=Requirement(
            database_fields=requires,
            dataframe_columns=requires
        )
    def inner_decorator(fn: FilterFunction) -> Filter:
        return Filter(fn=fn, name=fn.__name__, requires=requires, level=level)
    return inner_decorator

f = filter_fun
r = Requirement

@f(requires='sequence', level='individual')
def is_not_half_tryptic(full_sequence: str) -> bool:
    return full_sequence[0] in ['K','R','-'] and (full_sequence[-3] in ['K','R'] or full_sequence[-1] == '-')

@f(requires='sequence', level='individual')
def is_not_terminal(full_sequence: str) -> bool:
    return '-' not in full_sequence

@f(requires='uniprot', level='individual')
def is_not_reverse(uniprot: str) -> bool:
    return not uniprot.startswith('Reverse_')

@f(requires='clean_sequence', level='individual')
def is_not_miscleaved(clean_sequence: str) -> bool:
    seq = clean_sequence[:-1]
    return not('K' in seq or 'R' in seq)

@f(requires='clean_sequence', level='individual')
def is_not_miscleaved_more_than_once(clean_sequence: str) -> bool:
    seq = clean_sequence[:-1]
    return (seq.count('K') + seq.count('R')) < 2

@f(requires='unique_sequence', level='individual')
def is_unique_sequence(unique_sequence: bool) -> bool:
    return unique_sequence

@f(requires='unique_sequence', level='individual')
def is_non_unique_sequence(unique_sequence: bool) -> bool:
    return not unique_sequence

@f(requires='clean_sequence', level=['aggregate'])
def min_num_unique_sequences(clean_sequences: List[str], min_unique_sequences: int=2) -> bool:
    no_filter = [True] * len(clean_sequences)
    filter_all = [False] * len(clean_sequences)
    return no_filter if len(set(utils.flatten(clean_sequences))) >= min_unique_sequences else filter_all
