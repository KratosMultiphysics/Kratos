import warnings

def checkInitialisationMC(XMCAlgorithm):
    """
    Method checking all attributes of different classes are correctly set to run Monte Carlo algorithm (both standard and asynchronous).
    """

    solverWrapperDictionary = XMCAlgorithm.monteCarloSampler.indexConstructorDictionary[
            "samplerInputDictionary"]["solverWrapperInputDictionary"]
    positionMaxNumberIterationsCriterion = XMCAlgorithm.positionMaxNumberIterationsCriterion
    tolerances = XMCAlgorithm.stoppingCriterion.tolerances()
    # perform checks
    checkInitialisationSolverWrapper(solverWrapperDictionary)
    if ("asynchronous" in solverWrapperDictionary):
        if (solverWrapperDictionary["asynchronous"] is True):
            checkMaxNumberIterationsCriterion(positionMaxNumberIterationsCriterion,tolerances)

def checkInitialisationCMLMC():
    pass

def checkInitialisationAMLMC():
    pass

def checkInitialisationMLMC(XMCAlgorithm):
    """
    Method checking all attributes of different classes are correctly set to run Multilevel Monte Carlo algorithm (both standard and asynchronous).
    """

    solverWrapperDictionary = XMCAlgorithm.monteCarloSampler.indexConstructorDictionary[
            "samplerInputDictionary"]["solverWrapperInputDictionary"]
    positionMaxNumberIterationsCriterion=XMCAlgorithm.positionMaxNumberIterationsCriterion
    tolerances=XMCAlgorithm.stoppingCriterion.tolerances()
    # perform checks
    checkInitialisationSolverWrapper(solverWrapperDictionary)
    if ("asynchronous" in solverWrapperDictionary):
        if (solverWrapperDictionary["asynchronous"] is True):
            checkMaxNumberIterationsCriterion(positionMaxNumberIterationsCriterion,tolerances)

def checkInitialisationSolverWrapper(solverWrapperDictionary):
    """
    Method checking solver wrapper keys are correctly set. A required key is outputBatchSize, which is needed to organize the quantities of interest into sublists (of future objects). For this reason, outputBatchSize must be smaller than or equal to the total number of quantities of interest. It is required if areSamplesSplit() is true.
    """

    obs = 1 # default value assigned in the solver wrapper

    if "outputBatchSize" in solverWrapperDictionary:
        obs = solverWrapperDictionary["outputBatchSize"]
    else:
        warnings.warn("outputBatchSize not defined in solverWrapper dictionary. The default value of 1 will be considered.")

    if "qoiEstimator" in solverWrapperDictionary:
        numqoi = len(solverWrapperDictionary["qoiEstimator"])
        qoi_estimator = solverWrapperDictionary["qoiEstimator"]
    else:
        raise Exception("qoiEstimator not defined in solverWrapper dictionary. This entry is required by the SolverWrapper.")

    if "xmc.momentEstimator.MultiMomentEstimator" in qoi_estimator or "xmc.momentEstimator.MultiCombinedMomentEstimator" in qoi_estimator:
        if not "sizeMultiXMomentEstimator" in solverWrapperDictionary:
            raise Exception("solverWrapperDictionary does not contain the key sizeMultiXMomentEstimator. The use of MultiMomentEstimator and MultiCombinedMomentEstimator requires to specify the dimension of such multi-dimensional quantities of interest in the solverWrapperDictionary.")

    if obs > (numqoi):
        raise Exception ("solverWrapper: outputBatchSize exceeding maximum dimension. Set a value <= number of estimators, namely len(qoiEstimator).")

    if solverWrapperDictionary["refinementStrategy"] != "reading_from_file":
        if not "refinementParametersPath" in solverWrapperDictionary:
            raise Exception("solverWrapperDictionary does not contain the key \"refinementParametersPath\". This key is required for running with a \"refinementStrategy\" different from \"reading_from_file\".")

def checkMaxNumberIterationsCriterion(positionMaxNumberIterationsCriterion,tolerances):
    """
    Method checking maximum number of iterations is set and is the last entry of MonoCriteriaDictionary, which is required for consistency.
    """

    if (positionMaxNumberIterationsCriterion is not None):
        if (len(tolerances) == (positionMaxNumberIterationsCriterion + 1)):
            pass
        else:
            raise Exception ("Number of monoCriteria defined in monoCriteriaDictionary and positionMaxNumberIterationsCriterion defined in xmcAlgorithmDictionary not consistent. Should be positionMaxNumberIterationsCriterion = number monoCriteria and the monoCriterion defining the maximum number of iterations set as last entry of monoCriteriaDictionary.")
    else:
        raise Exception ("positionMaxNumberIterationsCriterion not set in xmcDictionary. Set it in order to run the asynchronous framework.")