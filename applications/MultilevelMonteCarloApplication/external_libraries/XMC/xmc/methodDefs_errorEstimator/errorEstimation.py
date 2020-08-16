import numpy as np
from xmc.tools import normalInverseCDF

# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

#TODO - think of potentially better name
def errorEstimationStatError(cdfValue, globalEstimations):
    """
    Accept the summation over the variances divided by number of samples over all
    levels and return its square root as the statistical error
    """
    if cdfValue is None:
        cdfValue = 1 # default behavior
    else:
        cdfValue = normalInverseCDF(cdfValue[0]) # [Pisaroni et al.,CMLMC,pag.25]
    error = cdfValue*np.sqrt(globalEstimations[0])
    # Ensure no NumPy type
    return float(error)

@ExaquteTask(returns=1, estimations={Type: COLLECTION_IN, Depth: 1})
def errorEstimationStatError_Task(parameters, estimations):
    """
    Same as errorEstimationStatError but is a pycompss task
    """
    return errorEstimationStatError(parameters, estimations)

def errorEstimationTerr(cdfValue,globalEstimations):
    """
    Return total error given moments/variances
    """
    error = globalEstimations[0] + errorEstimationStatError(cdfValue, globalEstimations[1:])
    # Ensure no NumPy type
    return float(error)

@ExaquteTask(returns=1, estimations={Type: COLLECTION_IN, Depth: 1})
def errorEstimationTerr_Task(parameters, estimations):
    """
    Same as errorEstimationTerr but is a pycompss task
    """
    return errorEstimationTerr(parameters, estimations)

# TODO This should be *squared* error, according to the name
def errorEstimationMSE(_, globalEstimations):
    """
    Return mean squared error given moments/variances
    """
    error = np.sqrt( globalEstimations[0]**2 + globalEstimations[1] )
    # Ensure no NumPy type
    return float(error)

@ExaquteTask(returns=1, estimations={Type: COLLECTION_IN, Depth: 1})
def errorEstimationMSE_Task(parameters, estimations):
    """
    Same as errorEstimationMSE but is a PyCOMPSs task
    """
    return errorEstimationMSE(parameters, estimations)
