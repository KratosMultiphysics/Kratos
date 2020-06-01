# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

import numpy as np
import xmc.solverWrapper as sw
import time

class MultiLevelRNGSolverWrapper(sw.SolverWrapper):
    """
    solveWrapper type whose solve method accepts a standard normal 
    random number and scales it by an index-dependent mean and 
    variance

    Constructor arguements - 
    solverWrapperIndex - Index-space position of the solverWrapper instance
    rates - list of length 4 whose first and last two entries are respectively
    the constant and rate required to compute the mean and variance using the 
    equation mean (or) variance = constant * exp (-rate*level)
    """
    # TODO update documentation above ('rates' is of length 6)
    def __init__(self, **keywordArgs):
        super().__init__(**keywordArgs)
        self.mean = self.parameters[0] * np.exp(-(self.parameters[1]*self.solverWrapperIndex[0]))
        self.standardDeviation = 1 + self.parameters[2] * np.exp(-(self.parameters[3]*self.solverWrapperIndex[0]))
        self.waitForTime = self.parameters[4] * np.exp(-(self.parameters[5]*self.solverWrapperIndex[0]))

    @ExaquteTask(returns=2)
    def solveForOneQoI_Task(self,randomInput):
        qoi = None
        if all([component>=0 for component in self.solverWrapperIndex]):
            qoi = self.mean + self.standardDeviation*randomInput[0]
        else:
            qoi = 0
        return qoi, self.waitForTime

    def solve(self, randomInput):
        # TODO Change hard coding
        number_of_qoi = 1
        # TODO Everything will work if the solve method returns a list of QoI-futures
        # and an associated time-future that it took to compute the list
        qoi,time_for_all_qoi = self.solveForOneQoI_Task(randomInput)
        list_of_qoi = [qoi]*number_of_qoi
        return list_of_qoi,time_for_all_qoi
