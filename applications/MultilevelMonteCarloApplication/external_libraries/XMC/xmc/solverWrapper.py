"""Module for interfaces with external solvers."""

from typing import Union, Tuple


class SolverWrapper:
    """Class used to call the relevant solver.

    Since this is very problem- and tool-specific, sub-classes (inheriting from that one) should be implemented for each solver. Therefore, this is a template of what SampleGenerator objects need.
    """

    # TODO should be renamed 'index'
    solverWrapperIndex: Tuple[int]
    """Place in the hierarchy of solvers."""
    outputDimension: Union[int, Tuple[int]]
    """Dimension of the first output of :py:meth:`.solve`.
    If integer, it is equals to its length.
    If list of integers, then it means that samples is split in sub-lists, and ``outputDimension`` is ``[len(subSample) for subSample in sample]``."""

    def __init__(self, **keywordArgs):
        """Constructor of SolverWrapper instances.

        :keyword tuple[int] index: Sets :py:attr:`~.solverWrapperIndex`.
        :keyword outputDimension: Sets :py:attr:`~.outputDimension`.
        :type outputDimension: int or list[int]
        """
        self.solverWrapperIndex = keywordArgs.get("index")
        self.outputDimension = keywordArgs.get("outputDimension")

    def solve(self, randomEvent) -> Tuple[list, float]:
        """
        Input: random event in the form expect by the solver wrapper.

        Outputs:
        1. sample: list of structure respecting self.outputDimension
        2. cost: float, measure of time to generate the sample.
        """
        pass
