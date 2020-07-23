# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

import xmc.tools

#TODO rename; e.g. 'sum'
def assembleValue(hierarchy, indexEstimations):
    """
    Compute a simple sum over all index-wise contributions
    """
    # Flatten nested list and sum
    return sum([est for estList in indexEstimations for est in estList])

@ExaquteTask(returns=1, indexEstimations={Type: COLLECTION_IN, Depth: 2})
def assembleValue_Task(hierarchy, indexEstimations):
    """
    Same as assembleValue, but is a PYCOMPSS task
    """
    return assembleValue(hierarchy, indexEstimations)

def assembleBias(hierarchy, indexEstimations):
    """
    Based on hierarchy, compute the sum over the index estimations that lie on
    the boundary of the index space as an estimate for the bias.
    Note - The procedure used in this function to find the boundary of
    a given hierarchy assumes that the index set is covex
    """
    # Select estimations of first random variable
    indexEstimations = indexEstimations[0]
    # Compute maxima along each dimension of the hierarchy
    index_set,_ = xmc.tools.splitOneListIntoTwo(hierarchy)
    bias_booleans = xmc.tools.strictlyPositiveBoundaryBooleans(index_set)
    bias_estimations = []
    for i in range(len(hierarchy)):
        if bias_booleans[i] is True:
            bias_estimations.append(indexEstimations[i])
    return abs(sum(bias_estimations))

@ExaquteTask(returns=1, indexEstimations={Type: COLLECTION_IN, Depth: 2})
def assembleBias_Task(hierarchy, indexEstimations):
    """
    Same as assembleBias, but is a PyCOMPSs task
    """
    return assembleBias(hierarchy, indexEstimations)

#TODO Duplicate of assembleValue
#TODO misnamed: variance not divived by number of samples
def assembleStatisticalError(hierarchy, indexEstimations):
    """
    Add together the variance over all indices
    to estimate global variance
    """
    return assembleValue(hierarchy, indexEstimations)

@ExaquteTask(returns=1, indexEstimations={Type: COLLECTION_IN, Depth: 2})
def assembleStatisticalError_Task(hierarchy, indexEstimations):
    """
    Same as assembleStatisticalError but is a PYCOMPSS task
    """
    return assembleStatisticalError(hierarchy, indexEstimations)

def assembleInterpolationError(hierarchy, indexEstimations):
    pass
