# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

# Import XMC methods
from xmc.tools import dynamicImport
from xmc.tools import instantiateObject
from xmc.tools import sum_Task
from xmc.tools import unpackedList

class SampleGenerator():
    """
    Class used for sample generation. Produces only one set of correlated samples at a time.
    """

    def __init__(self, **keywordArgs):
        # Attributes
        self.randomGenerator = instantiateObject(keywordArgs.get('randomGenerator'),**keywordArgs.get('randomGeneratorInputDict'))
        self.qoiProcessor = dynamicImport(keywordArgs.get('qoiProcessor', 'xmc.tools.returnInput'))

        solver_wrapper_indices = keywordArgs.get('solverWrapperIndices')
        self.solvers = []
        for sw_index in solver_wrapper_indices:
            keywordArgs['solverWrapperInputDict']['index'] = sw_index
            self.solvers.append(instantiateObject(keywordArgs.get('solverWrapper'),**keywordArgs.get('solverWrapperInputDict')))

    def qoiFromRaw(self, rawSolutions):
        """
        Function that generates additionally `QoI` from given QoI. For example,
        if one wished to study the statisics of the product or sum or difference
        of two or more QoI
        """
        return rawSolutions

    def randomInput(self):
        # TODO rewrite this method based on a more generic definition
        random_inputs = [self.randomGenerator.realisation()]*len(self.solvers)
        return random_inputs

    def rawSolutions(self, randomInputs):
        """
        Return the array of raw solutions produced by solvers.
        """
        raw_solutions = []
        times_for_qoi = []
        for i in range(len(self.solvers)):
            raw_solution,time_for_qoi = self.solvers[i].solve(randomInputs[i])
            raw_solutions.append(raw_solution)
            times_for_qoi.append(time_for_qoi)
        return raw_solutions, times_for_qoi

    def rearrangeListOfQoI(self,listOfQoI):
        return list(map(list, zip(*listOfQoI)))

    def generate(self, indexValue):
        """
        Generate new QoI instance. Calls randomInput (if no evaluation point is provided), then rawSolutions, then pass the outputs to qoiFromRaw and return its output.
        """
        # Generate the random input necessary for the solver
        random_inputs = self.randomInput()

        # Generate the list [[QoI_1^l,...,QoI_N^l],[QoI_1^(l-1),...,QoI_N^(l-1)],...]
        raw_solutions, times_for_qoi = self.rawSolutions(random_inputs)

        # Do additional processing
        list_of_qoi = self.qoiFromRaw(raw_solutions)

        # Rearrange list to generate [[QoI_1^l, QoI_1^(l-1), ...],...,[QoI_N^l, QoI_N^(l-1), ...]]
        list_of_qoi = self.rearrangeListOfQoI(list_of_qoi)

        # Add together times
        total_time = sum_Task(*unpackedList(times_for_qoi))

        return list_of_qoi, total_time

    def numberOfOutputs(self,solverChoice=0):
        #TODO Should be named numberOfOutputs instead of number_qoi
        return self.solvers[solverChoice].number_qoi

    def numberOfCombinedOutputs(self,solverChoice=0):
        return self.solvers[solverChoice].number_combined_qoi