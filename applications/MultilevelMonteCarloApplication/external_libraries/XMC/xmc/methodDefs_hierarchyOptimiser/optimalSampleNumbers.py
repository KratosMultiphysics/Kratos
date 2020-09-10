import numpy as np

def singleIndexDoubleSampleNumber(inputDict, newLevels):
    oldHierarchy = inputDict['oldHierarchy']
    old_sample_number = oldHierarchy[0][1]
    return [2*old_sample_number]

def singleIndexConstantSampleNumber(inputDict, newLevels):
    defaultHierarchy = inputDict['defaultHierarchy']
    sample_number = defaultHierarchy[0][1]
    return [sample_number]

def singleIndexAdaptive(inputDict, newLevels):
    oldHierarchy = inputDict['oldHierarchy']
    variance = inputDict['estimations'][0][0]*oldHierarchy[0][1]
    tolerance = inputDict['tolerances'][0]
    return [int(np.ceil(variance/tolerance)**2)]

def multiLevelDoubleAllSamples(inputDict, newLevels):
    """
    Returns a list of sample numbers of same length as the number of entries
    in newLevels. Doubles the number of samples from oldHierarchy if an
    entry of newLevels exists in oldHierarchy. If not, allocate a default
    newSampleNumber to the entry.
    """
    oldHierarchy = inputDict['oldHierarchy']
    newSampleNumber = inputDict['newSampleNumber']
    new_samples = []
    for newLevel in newLevels:
        is_level_found = False
        for oldElement in oldHierarchy:
            if newLevel==oldElement[0]:
                new_samples.append(2*oldElement[1])
                is_level_found = True
                break
        if is_level_found is False:
            new_samples.append(newSampleNumber)
    return new_samples

def multiLevelConstantSampleNumber(inputDict, newLevels):
    """
    Returns a list of sample numbers of same length as the number of levels of
    deault hierarchy. Keeps constant the number of samples from defaultHierarchy if an
    entry of newLevels exists in deafultHierarchy. If not, allocate a default
    newSampleNumber to the entry.
    """
    defaultHierarchy = inputDict['defaultHierarchy']
    newSampleNumber = inputDict['newSampleNumber']
    new_samples = []
    for level in newLevels:
        is_level_found = False
        for defaultElement in defaultHierarchy:
            if (level == defaultElement[0]):
                new_samples.append(1*defaultElement[1])
                is_level_found = True
                break
        if (is_level_found is False):
            new_samples.append(newSampleNumber)
    return new_samples

def multiLevelCostMinimisationTErr(inputDict, newLevels):
    oldHierarchy = inputDict['oldHierarchy']
    errorParameters = inputDict['errorParameters']
    cdf_value = errorParameters[0][0]
    tolerance = inputDict['tolerances'][0]
    splitting_parameter = inputDict['splittingParameter']
    cost_model = inputDict['costModel']
    cost_parameters = inputDict['costParameters']
    variances_for_hierarchy = inputDict['variancesForHierarchy']

    new_samples = []
    max_level = max([newLevels[i][0] for i in range(len(newLevels))])

    # Compute this term - sum_{k=0}^L sqrt(C_k V_k) for later use
    cost_variance_sum = 0.0
    for i in range(len(newLevels)):
        cost_level = cost_model(cost_parameters,newLevels[i])
        variance_level = variances_for_hierarchy[i]
        cost_variance_sum += np.sqrt(cost_level*variance_level)

    # Compute number of samples
    constantTerm = (cdf_value/(splitting_parameter*tolerance))**2 * cost_variance_sum
    for i in range(len(newLevels)):
        cost_level = cost_model(cost_parameters,newLevels[i])
        variance_level = variances_for_hierarchy[i]
        new_sample = constantTerm * np.sqrt(variance_level/cost_level)
        new_samples.append(int(np.ceil(new_sample)))
    return new_samples

def multiLevelCostMinimisationMSE(inputDict,newLevels):
    oldHierarchy = inputDict['oldHierarchy']
    errorParameters = inputDict['errorParameters']
    cdf_value = errorParameters[0][0]
    splitting_parameter = inputDict['splittingParameter']
    tolerance = inputDict['tolerances'][0]
    cost_model = inputDict['costModel']
    cost_parameters = inputDict['costParameters']
    cost_estimations = inputDict['costEstimations']
    variances_for_hierarchy = inputDict['variancesForHierarchy']

    new_samples = []
    max_level = max([newLevels[i][0] for i in range(len(newLevels))])

    # Compute this term - sum_{k=0}^L sqrt(C_k V_k) for later use
    cost_variance_sum = 0.0
    for i in range(len(newLevels)):
        cost_level = cost_model(cost_parameters,newLevels[i])
        variance_level = variances_for_hierarchy[i]
        cost_variance_sum += np.sqrt(cost_level*variance_level)

    # Compute number of samples
    for i in range(len(newLevels)):
        cost_level = cost_model(cost_parameters,newLevels[i])
        variance_level = variances_for_hierarchy[i]
        new_sample = (cdf_value/(splitting_parameter*tolerance))**2 * cost_variance_sum * np.sqrt(variance_level/cost_level)
        new_samples.append(int(np.ceil(new_sample)))
    return new_samples
