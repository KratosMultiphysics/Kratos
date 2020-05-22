# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3
import numpy as np
import time as time

import xmc.solverWrapper as sw

class SingleLevelRNGSolverWrapper(sw.SolverWrapper):
    """
    solveWrapper type whose solve method accepts a random number
    and returns it in a list of size the number of QoI

    Constructor arguements - 
    solverWrapperIndex - Index-space position of the solverWrapper instance
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # TODO - this is a temporary solution for changing futures list
    # to list of futures
    @ExaquteTask(returns=2)
    def solveForOneQoI_Task(self,randomInput):
        qoi = None
        start_time = time.time()
        if all([component>=0 for component in self.solverWrapperIndex]):
            qoi = randomInput[0]
        else:
            qoi = 0
        end_time = time.time()
        return qoi, (end_time-start_time)

    def solve(self, randomInput):
        # TODO Change hard coding
        number_of_qoi = 1
        # TODO Everything will work if the solve method returns a list of QoI-futures
        # and an associated time-future that it took to compute the list
        qoi,time_for_all_qoi = self.solveForOneQoI_Task(randomInput)
        list_of_qoi = [qoi]*number_of_qoi
        return list_of_qoi,time_for_all_qoi
