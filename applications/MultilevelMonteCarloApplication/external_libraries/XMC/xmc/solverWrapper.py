"""Module for interfaces with external solvers"""


class SolverWrapper:
    """
    Class used to call the relevant solver. Since this is very problem- and tool-
    specific, sub-classes (inheriting from that one) should be implement for each
    solver. Therefore, this is a template of what SampleGenerator objects need.

    Attributes:
    - solverWrapperIndex: list, place in the index set. Useful for the constructor.
    - outputDimension: integer or list of integer. If integer, equals to len(sample),
    where sample is the first output argument of self.solve(). If list of integers,
    then it means that samples is split in future lists, and outputDimension is
    [len(subSample) for subSample in sample].

    Methods:
    - solve. See its documentation.
    """

    def __init__(self, **keywordArgs):
        # TODO should be renamed 'index'
        self.solverWrapperIndex = keywordArgs.get("index")
        self.outputDimension = keywordArgs.get("outputDimension")

    def solve(self, randomEvent):
        """
        Input: random event in the form expect by the solver wrapper.

        Outputs:
        1. sample: list of structure respecting self.outputDimension
        2. cost: float, measure of time to generate the sample.
        """
        pass
