from .shared import *
from .shared import filter_fun as f, Requirement as r

from typing import List
from copy import deepcopy
import statistics

RatioList = List[float]

@f(requires='ratio', level='individual')
def is_quantified(ratio: float) -> bool:
    return ratio > 0

@f(requires=r(['ratio', 'num_ms2'], ['ratio', 'num_ms2']), level='individual')
def is_not_singleton(ratio: float, num_ms2: int) -> bool:
    return ratio != 20 or num_ms2 > 1

@f(requires='rsquared', level='individual')
def passes_rsquared_cutoff(rsquared: float, min_rsquared: float = 0.5):
    return rsquared >= min_rsquared

@f(requires='ratio', level=['experiment', 'aggregate'])
def filter_20s(ratios: RatioList, ratio_cutoff: int = 2) -> List[bool]:
    """Filter out erroneous seeming 20s."""
    # 20s are stripped if the following conditions are met:
    # - the set of ratios is not just composed of 0s and 20s
    # - there is only one 20
    # - the lowest non-zero, non-20 value is below a cutoff
    no_filter = [True] * len(ratios)
    if 20 not in ratios or ratios.count(20) > 1:
        return no_filter

    only_zeros_and_20s = len([x for x in ratios if x > 0 and x != 20]) == 0
    min_below_cutoff = min(ratios) < ratio_cutoff

    if not only_zeros_and_20s and min_below_cutoff:
        return [r != 20 for r in ratios]
    else:
        return no_filter


@f(requires='ratio', level=['experiment', 'aggregate'])
def filter_two_or_more_20s(ratios: RatioList, ratio_cutoff=2) -> List[bool]:
    """Filter out erroneous seeming 20s when there is more than one 20 value."""
    # 20s are stripped if the following conditions are met:
    # - the set of ratios is not just composed of 0s and 20s
    # - there are at least 2 20 ratios
    # - the lowest non-zero, non-20 value is below a cutoff
    no_filter = [True] * len(ratios)
    if 20 not in ratios or ratios.count(20) < 2:
        return no_filter

    only_zeros_and_20s = len([x for x in ratios if x > 0 and x != 20]) == 0
    min_below_cutoff = min(ratios) < ratio_cutoff

    if not only_zeros_and_20s and min_below_cutoff:
        return [r != 20 for r in ratios]
    else:
        return no_filter

@f(requires='ratio', level=['experiment', 'aggregate'])
def filter_by_stdev(ratios: RatioList, stdev_cutoff: float = 0.6, ratio_cutoff: float = 4.0) -> List[bool]:
    """Filter ratios based on a standard deviation cutoff."""
    no_filter = [True] * len(ratios)
    quantified_ratios = [r for r in ratios if r > 0]

    # no-op if less than two ratios
    if len(quantified_ratios) < 2:
        return no_filter

    stdev = statistics.stdev(quantified_ratios)
    mean = statistics.mean(quantified_ratios)
    minimum = min(quantified_ratios)

    # if stdev is tight enough relative to the mean
    # or if all the ratios are above a certain threshold
    # return unchanged
    # else, set everything but the minimum to zero
    if (stdev != None and stdev / mean < stdev_cutoff) or minimum > ratio_cutoff:
        return no_filter
    else:
        return [r == minimum for r in ratios]

@f(requires='ratio', level=['experiment', 'aggregate'])
def filter_20s_by_stdev(ratios: RatioList, stdev_cutoff: float = 0.6, ratio_cutoff: float = 4) -> List[bool]:
    """Filter out 20s from sets of ratios with high standard deviation."""
    no_filter = [True] * len(ratios)
    quantified_ratios = [r for r in ratios if r > 0]

    # require that there are at least two ratios
    # that there is a 20 in the set of ratios
    # and that the minimum ratio is below the ratio_cutoff
    # otherwise do nothing
    if len(quantified_ratios) < 2 or 20 not in quantified_ratios or min(quantified_ratios) > ratio_cutoff:
        return no_filter

    above_threshold = statistics.stdev(quantified_ratios) / statistics.mean(quantified_ratios) > stdev_cutoff

    if above_threshold:
        return [r != 20 for r in ratios]
    else:
        return no_filter


@f(requires='ratio', level='aggregate')
def min_num_datasets(ratios: RatioList, min_datasets: int = 2) -> List[bool]:
    no_filter = [True] * len(ratios)
    filter_all = [False] * len(ratios)
    quantified_ratios = len([x for x in ratios if x > 0])
    return no_filter if quantified_ratios >= min_datasets else filter_all

@f(requires='clean_sequence', level=['experiment', 'aggregate'])
def min_num_peptides(clean_sequences: List[str], min_peptides: int=2) -> bool:
    no_filter = [True] * len(clean_sequences)
    filter_all = [False] * len(clean_sequences)
    return no_filter if len(list(utils.flatten(clean_sequences))) >= min_peptides else filter_all
