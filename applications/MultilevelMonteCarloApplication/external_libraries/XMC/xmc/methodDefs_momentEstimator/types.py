"""Provide types used in the xmc.momentEstimator module"""

# Define types
from typing import List, DefaultDict, Union, Iterable

PowerSumsDict = DefaultDict[str, Union[float, List[float]]]
SampleArray = List[List[List[float]]]
HStatistics = Union[float, List[float]]
ListIndex = Union[int, slice, Iterable[int]]
PowerSumsDictUL = DefaultDict[str, PowerSumsDict]  # UL: upper, lower
CombinedSampleArray = List[List[Union[List[float], int]]]
# TODO Account for NumPy types
