def noBlending(parameters,newLevels,inputDict):
    """
    Return list of variances based purely on the variance model
    """
    old_hierarchy = inputDict['oldHierarchy']
    variance_model = inputDict['models'][1]
    variance_parameters = inputDict['parametersForModel'][1]

    variance_list = []
    for level in newLevels:
        variance_list.append( variance_model(variance_parameters,level) )
    return variance_list

def bayesianUpdate(parameters, newLevels, inputDict):
    """
    Use bayesian estimation to blend the variance estimations and variance
    model such that the value of the estimation is used if there are a large
    number of samples and the value of the model is used if there is a small
    number of samples.
    """
    k0 = parameters[0]
    k1 = parameters[1]

    old_hierarchy = inputDict['oldHierarchy']
    bias_model = inputDict['models'][0]
    bias_parameters = inputDict['parametersForModel'][0]
    bias_estimations = inputDict['estimations'][0]
    variance_model = inputDict['models'][1]
    variance_parameters = inputDict['parametersForModel'][1]
    variance_estimations = [inputDict['estimations'][1][i]*old_hierarchy[i][1]
                            for i in range(len(old_hierarchy))]

    bayesianVariance = []
    for level in newLevels:
        number_samples = 0
        level_in_old_hierarchy = 0
        for position in range(len(old_hierarchy)):
            element = old_hierarchy[position]
            if(element[0]==level):
                number_samples = element[1]
                break
        model_bias = bias_model(bias_parameters,level)
        model_variance = variance_model(variance_parameters,level)

        Xi1 = 0.5 + k1/model_variance + number_samples * 0.5
        Xi2 = k1 + (number_samples-1)*0.5*variance_estimations[position] +\
        0.5 * k0 * number_samples * (bias_estimations[position]-model_bias)**2 /\
        (k0 + number_samples)

        bayesianVariance.append(Xi2/(Xi1-0.5))
    return bayesianVariance
