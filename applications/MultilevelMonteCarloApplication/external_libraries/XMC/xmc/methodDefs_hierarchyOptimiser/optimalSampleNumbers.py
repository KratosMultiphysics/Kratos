import numpy as np
from xmc.tools import normalInverseCDF


def singleIndexDoubleSampleNumber(inputDict, newLevels):
    oldHierarchy = inputDict["oldHierarchy"]
    old_sample_number = oldHierarchy[0][1]
    return [2 * old_sample_number]


def singleIndexConstantSampleNumber(inputDict, newLevels):
    defaultHierarchy = inputDict["defaultHierarchy"]
    sample_number = defaultHierarchy[0][1]
    return [sample_number]


def singleIndexAdaptive(inputDict, newLevels):
    oldHierarchy = inputDict["oldHierarchy"]
    variance = inputDict["estimations"][0][0] * oldHierarchy[0][1]
    tolerance = inputDict["tolerances"][0]
    # TODO Here the tolerance is applied to the standard deviation. Check consistency with the rest of the library.
    return [int(np.ceil(variance / tolerance ** 2))]


def multiLevelDoubleAllSamples(inputDict, newLevels):
    """
    Returns a list of sample numbers of same length as the number of entries
    in newLevels. Doubles the number of samples from oldHierarchy if an
    entry of newLevels exists in oldHierarchy. If not, allocate a default
    newSampleNumber to the entry.
    """
    oldHierarchy = inputDict["oldHierarchy"]
    newSampleNumber = inputDict["newSampleNumber"]
    new_samples = []
    for newLevel in newLevels:
        is_level_found = False
        for oldElement in oldHierarchy:
            if newLevel == oldElement[0]:
                new_samples.append(2 * oldElement[1])
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
    defaultHierarchy = inputDict["defaultHierarchy"]
    newSampleNumber = inputDict["newSampleNumber"]
    new_samples = []
    for level in newLevels:
        is_level_found = False
        for defaultElement in defaultHierarchy:
            if level == defaultElement[0]:
                new_samples.append(1 * defaultElement[1])
                is_level_found = True
                break
        if is_level_found is False:
            new_samples.append(newSampleNumber)
    return new_samples


def minimalCostStatisticalError(inputDict, newLevels):
    """Compute sample numbers to satisfy tolerance on statistical error with minimal cost."""
    # Variances
    variances = inputDict["blendedVariances"]
    if variances is None:
        # Variance was not blended
        hierarchy = inputDict["oldHierarchy"]
        if inputDict["parametersForModel"] is not None:
            # We have a valid model
            variance_params = inputDict["parametersForModel"][1]
            variance_model = inputDict["models"][1]
            variances = [variance_model(variance_params, l) for l in newLevels]
            # Use sample estimation for first level
            variances[0] = inputDict["estimations"][-1][0] * hierarchy[0][1]
        else:
            # No model: fall back to sample estimations
            variances = [v * n[1] for v, n in zip(inputDict["estimations"][-1], hierarchy)]

    # Costs
    cost_parameters = inputDict["costParameters"]
    if cost_parameters is not None:
        # We have a valid model
        cost_model = inputDict["costModel"]
        costs = [cost_model(cost_parameters, l) for l in newLevels]
        # Use sample estimation for first level
        costs[0] = inputDict["costEstimations"][0]
    else:
        # No model: fall back to sample estimations
        costs = inputDict["costEstimations"]

    # Compute this term - sum_{k=0}^L sqrt(C_k V_k) for later use
    cost_variance_sum = 0.0
    for i, cost in enumerate(costs):
        cost_variance_sum += np.sqrt(cost * variances[i])

    # Factor for sample numbers, independent of level
    cdf_value = normalInverseCDF(inputDict["errorParameters"][0][0])
    tolerance = inputDict["tolerances"][0] * inputDict["splittingParameter"]
    constantFactor = (cdf_value / tolerance) ** 2 * cost_variance_sum

    # Compute number of samples
    new_samples = []
    for i in range(len(newLevels)):
        if i > len(costs) - 1:
            # No data for this level. Set sample number to zero.
            # The hierarchy optimiser enforces the minimal number of samples per level
            new_samples.append(0)
            continue
        variance_level = variances[i]
        new_sample_nb = constantFactor * np.sqrt(variance_level / costs[i])
        new_samples.append(int(np.ceil(new_sample_nb)))
    return new_samples
