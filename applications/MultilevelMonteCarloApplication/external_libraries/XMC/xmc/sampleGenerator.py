# Import XMC methods
from xmc.tools import dynamicImport, instantiateObject
from xmc.methodDefs_sampleGenerator import qoiProcessor


class SampleGenerator:
    """
    Class used for sample generation. Interacts with a sample generator and solvers to generate a set of correlated samples.
    """

    def __init__(self, **keywordArgs):
        # Attributes
        self.randomGenerator = instantiateObject(
            keywordArgs.get("randomGenerator"),
            **keywordArgs.get("randomGeneratorInputDictionary"),
        )
        self.qoiProcessor = dynamicImport(
            keywordArgs.get("qoiProcessor", "xmc.tools.returnInput")
        )

        solver_wrapper_indices = keywordArgs.get("solverWrapperIndices")
        self.solvers = []
        for sw_index in solver_wrapper_indices:
            keywordArgs["solverWrapperInputDictionary"]["index"] = sw_index
            self.solvers.append(
                instantiateObject(
                    keywordArgs.get("solverWrapper"),
                    **keywordArgs.get("solverWrapperInputDictionary"),
                )
            )

    def qoiFromRaw(self, rawSolutions):
        """
        Function that generates additionally `QoI` from given QoI. For example,
        if one wished to study the statisics of the product or sum or difference
        of two or more QoI
        """

        return rawSolutions

    def randomEvent(self):
        # TODO rewrite this method based on a more generic definition
        return self.randomGenerator.realisation()

    def rawSolutions(self, randomEvent):
        """
        Return the array of raw solutions and computation time produced by the solvers.

        Input argument: random event, as expect by the solve method of the solvers.

        Output argument:
        - solutions: list of solver outputs; solutions[j] is the solution returned by self.solvers[j].
        - times: list of computational times; times[j] is the time returned by self.solvers[j].

        See the documentation of SampleGenerator.generate for more details on the different possible
        data structures of solutions.
        """

        solutions = []
        times = []
        for solver in self.solvers:
            solution, resolution_time = solver.solve(randomEvent)
            solutions.append(solution)
            times.append(resolution_time)

        return solutions, times

    def generate(self, randomEvent=None):
        """
        Generate new single sample.
        Calls randomEvent (if no evaluation point is provided),
        then rawSolutions, then pass the outputs to qoiFromRaw and return its output.

        Input arguments:
        - randomEvent (optional): random event to be passed to solvers.
        By default, the solvers receive a random event generated internally by the method SampleGenerator.randomEvent.

        Output arguments:
        - sample: list of length the number of solvers. Each element of that list is either a list of length SampleGenerator.sampleDimension() or a list of future lists of length len(SampleGenerator.sampleSplitSizes()). For the length of these future lists, see SampleGenerator.sampleSplitSizes.
        - time: list of floats of length the number of solvers. Measurement of the cost of each solver call.

        Note on data structure of output sample
        ---------------------------------------

        Let us note S = sample, S_j = sample[j], etc.. All such objects are lists.
        S is a nested list of length len(SampleGenerator.solvers) and depth at least 2.
        Elements of S_j may be future objects.

        If solver outputs are not split, (i.e. not SampleGenerator.areSamplesSplit()), S_j is of length SampleGenerator.sampleDimension().
        Example. If len(SampleGenerator.solvers) = 2 and SampleGenerator.sampleDimension() = 3:
        S_1 = [s1_o1, s1_o2, s1_o3]   (s: solver, o: output); idem for S_2.

        If solver outputs are split, (i.e. SampleGenerator.areSamplesSplit()), S_j is of length len(SampleGenerator.sampleSplitSizes())
        and S_jk is of length SampleGenerator.sampleSplitSizes()[k].
        Example. If len(SampleGenerator.solvers) = 2 and SampleGenerator.sampleSplitSizes() = [3, 2]:
        S_1 = [ S_(1,1), S_(1,2) ] = [ [s1_o1, s1_o2, s1_o3], [s1_o4, s1_o5] ]
        (s: solver, o: output); idem for solver S_2.
        """

        # Generate the random input necessary for the solver
        if randomEvent is None:
            randomEvent = self.randomEvent()

        # Generate solutions
        raw_outputs, time = self.rawSolutions(randomEvent)
        # Do additional processing
        sample = self.qoiFromRaw(raw_outputs)
        return sample, time

    def sampleDimension(self, *args) -> int:
        """
        Returns the number of outputs of a single solver, as integer.
        The time measurement is not counted.
        This gives no information on whether the samples are split into
        future (sub-lists). Use SampleGenerator.areSamplesSplit for this.
        """

        sampleDim = self._solverOutputDimension(*args)
        # Ensure that sampleDim is iterable
        if isinstance(sampleDim, int):
            sampleDim = [sampleDim]
        # Returns the sum of all dimensions
        return sum(sampleDim)

    def sampleSplitSizes(self, *args):
        """
        Returns the size of lists in which solver outputs are split (None if they are not).
        Denoting s the first output argument of self.solvers[any].solve, splitSizes[j]=len(s[j]).

        Input arguments: same as SampleGenerator._solverOutputDimensions.

        Output arguments:
        - splitSizes: list of integers if SampleGenerator.areSamplesSplit(); None otherwise.
        """

        if self.areSamplesSplit(*args):
            return self._solverOutputDimension(*args)
        else:
            return None

    def areSamplesSplit(self, *args) -> bool:
        """
        Returns True if samples are split into future (sub-)lists, False otherwise.
        Accepts the same input arguments as _solverOutputDimensions.
        """

        sampleDim = self._solverOutputDimension(*args)
        return isinstance(sampleDim, (list, tuple))

    def _solverOutputDimension(self, solverChoice: int = 0):
        """
        Returns the dimension of the output of a solver.
        See SolverWrapper documentation for details.

        Input arguments:
        - solverChoice (optional): index of solver from which to fetch.
        Default is 0.
        """

        return self.solvers[solverChoice].outputDimension
