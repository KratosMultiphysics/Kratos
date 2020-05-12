import numpy as np
from xmc.tools import packedList
from xmc.tools import normalInverseCDF

# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

#TODO - think of potentially better name
def errorEstimationStatError(ignore, *args):
    """
    Accept the summation over the variances divided by number of samples over all
    levels and return its square root as the statistical error
    """
    assembledEstimationsList = packedList(args)
    # TODO - Think of better place for this assertion
    assert(len(assembledEstimationsList)==1),"length of assembledEstimationsList passed to errorEstimationStatError is not 1"
    if ignore is None:
        cdf_value = 1 # default behavior
    else:
        cdf_value = normalInverseCDF(ignore[0]) # [Pisaroni et al.,CMLMC,pag.25]
    return cdf_value*np.sqrt(assembledEstimationsList[0])

@ExaquteTask(returns=1)
def errorEstimationStatError_Task(*args):
    """
    Same as errorEstimationStatError but is a pycompss task
    """
    return errorEstimationStatError(*args)

def errorEstimationTerr(cdf_value,*args):
    """
    Return total error given moments/variances
    """
    assembledEstimationsList = packedList(args)
    return assembledEstimationsList[0] + normalInverseCDF(cdf_value[0]) * np.sqrt(assembledEstimationsList[1])

@ExaquteTask(returns=1)
def errorEstimationTerr_Task(*args):
    """
    Same as errorEstimationTerr but is a pycompss task
    """
    return errorEstimationTerr(*args)

def errorEstimationMSE(ignore, *args):
    """
    Return mean squared error given moments/variances
    """
    assembledEstimationsList = packedList(args)
    return np.sqrt( assembledEstimationsList[0]**2 + assembledEstimationsList[1] )

@ExaquteTask(returns=1)
def errorEstimationMSE_Task(*args):
    """
    Same as errorEstimationMSE but is a PyCOMPSs task
    """
    return errorEstimationMSE(*args)