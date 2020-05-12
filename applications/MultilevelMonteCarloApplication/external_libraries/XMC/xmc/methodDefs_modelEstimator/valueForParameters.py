import numpy as np

# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

def geometricModel(parameters,point):
    """
    Compute y = C*exp(-w.x) where C = parameters[0] and w = parameters[1:]
    """
    return parameters[0]*np.exp(sum([-parameters[i+1]*point[i] for i in range(len(point))]))

@ExaquteTask(returns=1)
def geometricModel_Task(parameters,point):
    """
    Compute y = C*exp(-w.x) where C = parameters[0] and w = parameters[1:]
    """
    return geometricModel(parameters,point)
