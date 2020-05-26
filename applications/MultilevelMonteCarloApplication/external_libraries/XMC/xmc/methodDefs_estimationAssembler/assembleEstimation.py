# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

import xmc.tools

# TODO - potentially improved name for this method
def assembleValue(hierarchy,*indexEstimationsList):
    """
    Compute a simple sum over all index-wise contributions
    """
    indexEstimationsList = xmc.tools.packedList(indexEstimationsList)
    return sum(indexEstimationsList)

@ExaquteTask(returns=1)
def assembleValue_Task(hierarchy,*indexEstimationsList):
    """
    Same as assembleValue, but is a PYCOMPSS task
    """
    return assembleValue(hierarchy,*indexEstimationsList)

def assembleBias(hierarchy,*indexEstimationsList):
    """
    Based on hierarchy, compute the sum over the index estimations that lie on
    the boundary of the index space as an estimate for the bias.
    Note - The procedure used in this function to find the boundary of
    a given hierarchy assumes that the index set is covex
    """
    indexEstimationsList = xmc.tools.packedList(indexEstimationsList)
    # Compute maxima along each dimension of the hierarchy
    index_set,_ = xmc.tools.splitOneListIntoTwo(hierarchy)
    bias_booleans = xmc.tools.strictlyPositiveBoundaryBooleans(index_set)
    bias_estimations = []
    for i in range(len(hierarchy)):
        if bias_booleans[i] is True:
            bias_estimations.append(indexEstimationsList[i])
    return abs(sum(bias_estimations))

@ExaquteTask(returns=1)
def assembleBias_Task(hierarchy,*indexEstimationsList):
    """
    Same as assembleBias, but is a PYCOMPSS task
    """
    return assembleBias(hierarchy,*indexEstimationsList)

def assembleStatisticalError(hierarchy,*indexEstimationsList):
    """
    Add together the variance divided by number of samples over all indices
    to estimate statistical error
    """
    indexEstimationsList = xmc.tools.packedList(indexEstimationsList)
    # TODO Quentin: I don't understand that check below, and it throws an error
    # even when indexEstimationsList is as expected. Is it a mistake? 
    # assert(len(indexEstimationsList)==1),"length of indexEstimationsList is not equal to 1"
    return sum( indexEstimationsList )

@ExaquteTask(returns=1)
def assembleStatisticalError_Task(hierarchy,*indexEstimationsList):
    """
    Same as assembleStatisticalError but is a PYCOMPSS task
    """
    return assembleStatisticalError(hierarchy,*indexEstimationsList)

def assembleInterpolationError(hierarchy,*indexEstimationsList):
    pass
