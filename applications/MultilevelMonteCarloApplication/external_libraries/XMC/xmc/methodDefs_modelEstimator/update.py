from xmc.distributedEnvironmentFramework import *

import numpy as np

def updatePredictorGeometric(data):
    """
    Fit a model of the form y = C*exp(-w.x) where x is the antecedent
    and y is the corresponding image.
    """
    points = [data[i][0] for i in range(len(data))]
    values = [data[i][1] for i in range(len(data))]

    # Setting up matrix A and vector b for least squares solution
    # of the problem min_x ||Ax-b||_2^2
    A = np.zeros((len(points), len(points[0])+1), dtype=float)
    b = np.zeros((len(points), 1), dtype=float)
    A[:,0]=1.

    for i in range(len(points)):
        A[i,1:] = points[i]
        b[i] = np.log(values[i])

    [x,_,_,_] = np.linalg.lstsq(A,b, rcond=-1)
    x = x.tolist()
    x = [x[i][0] for i in range(len(x))]
    parameters = []
    parameters.append(np.exp(x[0]))
    for i in (range(len(x)-1)):
        parameters.append(-x[i+1])
    return parameters

@ExaquteTask(returns=1, data={Type: COLLECTION_IN, Depth: 2})
def updatePredictorGeometric_Task(data):
    """
    Fit a model of the form y = C*exp(-w.x) where x is the antecedent
    and y is the corresponding image.
    """
    # TODO do we really need list packing and unpacking here?
    return updatePredictorGeometric(data)
