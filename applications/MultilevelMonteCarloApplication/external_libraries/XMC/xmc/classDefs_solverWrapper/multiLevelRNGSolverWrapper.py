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
        if self.outputDimension is None:
            self.outputDimension = 1
        self.parameters = keywordArgs.get('parameters')
        self.mean = self.parameters[0] * np.exp(-(self.parameters[1]*self.solverWrapperIndex[0]))
        self.standardDeviation = 1 + self.parameters[2] * np.exp(-(self.parameters[3]*self.solverWrapperIndex[0]))
        self.waitForTime = self.parameters[4] * np.exp(-(self.parameters[5]*self.solverWrapperIndex[0]))

    @ExaquteTask(returns=2)
    def _drawSample_Task(self,randomInput):
        sample = None
        if all([component>=0 for component in self.solverWrapperIndex]):
            sample = self.mean + self.standardDeviation*randomInput
        else:
            sample = 0
        return sample, self.waitForTime

    def solve(self, randomInput):
        # TODO is it necessary to call a task here?
        sample,totalTime = self._drawSample_Task(randomInput)
        # Make it a list of expected size
        if isinstance(self.outputDimension,int):
            sample = [sample]*self.outputDimension
        else: # list of integers
            sample = [[sample]*dim for dim in self.outputDimension] 
        return sample,totalTime
