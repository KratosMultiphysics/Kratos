class SolverWrapper():
    """
    Class used to call the relevant solver. Since this is very problem- and tool-
    specific, sub-classes (inheriting from that one) should be implement for each
    solver. Therefore, this is but a loose template of what SampleGenerator
    objects need.
    """

    def __init__(self, **keywordArgs):
        self.solverWrapperIndex = keywordArgs.get('index')
        self.outputDimension = keywordArgs.get('outputDimension',None)
        self.outputBatchSize = keywordArgs.get('outputBatchSize',None)
        self.parameters = keywordArgs.get('parameters',None)

    # There should be no definitions below this point.
    # SolverWrapper as a class is mostly empty and is inherited
    # by solver-specific-wrapper classes .eg. SolverWrapperKratos
    # or SolverWrapperOpenFOAM or SolverWrapperFluent

    def solve(self):
        """
        Outputs:
        1. List of length self.outputDimension (sample)
        2. Float (sample cost)
        """
        pass
