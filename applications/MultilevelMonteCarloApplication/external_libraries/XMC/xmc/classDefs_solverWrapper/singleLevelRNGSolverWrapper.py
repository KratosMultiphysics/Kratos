from xmc.distributedEnvironmentFramework import *

import numpy as np
import time as time

import xmc.solverWrapper as sw

class SingleLevelRNGSolverWrapper(sw.SolverWrapper):
    """
    solveWrapper type whose solve method accepts a random number
    and returns it in a list of size self.ouptutDimension

    Constructor arguments -
    solverWrapperIndex - Index-space position of the solverWrapper instance
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.outputDimension is None:
            self.outputDimension = 1

    # TODO - this is a temporary solution for changing futures list
    # to list of futures
    @ExaquteTask(returns=2)
    def _drawSample_Task(self,randomInput):
        start_time = time.time()
        if all([component>=0 for component in self.solverWrapperIndex]):
            sample = randomInput
        else:
            sample = 0
        end_time = time.time()
        return sample, (end_time-start_time)

    def solve(self, randomInput):
        # TODO is it necessary to call a task here?
        sample,totalTime = self._drawSample_Task(randomInput)
        sample = [sample]*self.outputDimension
        return sample,totalTime
