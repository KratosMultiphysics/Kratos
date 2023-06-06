import numpy as np

from exaqute import *

def geometricModel(parameters,point):
    """
    Compute y = C*exp(-w.x) where C = parameters[0] and w = parameters[1:]
    """
    return parameters[0]*np.exp(sum([-parameters[i+1]*point[i] for i in range(len(point))]))

@task(keep=True, returns=1)
def geometricModel_Task(parameters,point):
    """
    Compute y = C*exp(-w.x) where C = parameters[0] and w = parameters[1:]
    """
    return geometricModel(parameters,point)
