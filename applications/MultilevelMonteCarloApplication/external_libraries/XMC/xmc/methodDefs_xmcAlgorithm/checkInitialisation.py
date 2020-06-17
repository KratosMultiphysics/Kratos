def checkInitialisationMC(solverWrapperDictionary,positionMaxNumberIterationsCriterion=None,tolerances=None):
    checkInitialisationSolverWrapper(solverWrapperDictionary)
    if ("asynchronous" in solverWrapperDictionary):
        if (solverWrapperDictionary["asynchronous"] is True):
            checkMaxNumberIterationsCriterion(positionMaxNumberIterationsCriterion,tolerances)

def checkInitialisationCMLMC():
    pass

def checkInitialisationAMLMC():
    pass

def checkInitialisationMLMC(solverWrapperDictionary,positionMaxNumberIterationsCriterion=None,tolerances=None):
    checkInitialisationSolverWrapper(solverWrapperDictionary)
    if ("asynchronous" in solverWrapperDictionary):
        if (solverWrapperDictionary["asynchronous"] is True):
            checkMaxNumberIterationsCriterion(positionMaxNumberIterationsCriterion,tolerances)

def checkInitialisationSolverWrapper(solverWrapperDictionary):
    sql = 1 ; ncq = 0 ; nq = 1 # default values
    if "outputBatchSize" in solverWrapperDictionary:
        sql = solverWrapperDictionary["outputBatchSize"]
    if "numberCombinedQoi" in solverWrapperDictionary:
        ncq = solverWrapperDictionary["numberCombinedQoi"]
    if "numberQoI" in solverWrapperDictionary:
        nq = solverWrapperDictionary["numberQoI"]

    if (sql > (nq+ncq)):
        raise Exception ("solverWrapperDictionary: outputBatchSize exceeding maximum dimension. Set a value <= numberQoI + numberCombinedQoI.")

def checkMaxNumberIterationsCriterion(positionMaxNumberIterationsCriterion,tolerances):
    if (positionMaxNumberIterationsCriterion is not None):
        if (len(tolerances) == (positionMaxNumberIterationsCriterion + 1)):
            pass
        else:
            raise Exception ("Number of monoCriteria defined in monoCriteriaDictionary and positionMaxNumberIterationsCriterion defined in xmcAlgorithmDictionary not consistent. Should be positionMaxNumberIterationsCriterion = number monoCriteria and the monoCriterion defining the maximum number of iterations set as last entry of monoCriteriaDictionary.")
    else:
        raise Exception ("positionMaxNumberIterationsCriterion not set in xmcDictionary. Set it in order to run the asynchronous framework.")