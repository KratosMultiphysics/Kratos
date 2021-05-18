import numpy as np


def zeroDimensionSamplesOnly(inputDict):
    return [[]]

def maximalLevelsBiasModel(inputDict):
    oldHierarchy = inputDict['oldHierarchy']
    tolerance = inputDict['tolerances'][0]
    splitting_parameter = inputDict['splittingParameter']
    bias_parameters = inputDict['parametersForModel'][0]
    # Compute level for which bias model predicts
    optimal_level = np.ceil(-np.log((1-splitting_parameter)*tolerance/bias_parameters[0])/bias_parameters[1])

    # Extend hierarchy until that level is reached, or if optimal_level < max_level,
    # retain the same hierarchy
    new_levels = [oldHierarchy[i][0] for i in range(len(oldHierarchy))]
    max_level = max([new_levels[i][0] for i in range(len(new_levels))])
    while(max_level < optimal_level):
        new_levels.append([max_level+1])
        # recompute maximum
        max_level = max([new_levels[i][0] for i in range(len(new_levels))])
    return new_levels

def incrementLevelsByOne(inputDict):
    oldHierarchy = inputDict['oldHierarchy']
    new_levels = [oldHierarchy[i][0] for i in range(len(oldHierarchy))]
    max_level = max([new_levels[i][0] for i in range(len(new_levels))])
    new_levels.append([max_level+1])
    return new_levels

def constantNumberLevels(inputDict):
    defaultHierarchy = inputDict['defaultHierarchy']
    number_levels = [defaultHierarchy[i][0] for i in range(len(defaultHierarchy))]
    return number_levels

def fullIndexSet(inputDict):
    pass

def optimalIndexSet(inputDict):
    pass
