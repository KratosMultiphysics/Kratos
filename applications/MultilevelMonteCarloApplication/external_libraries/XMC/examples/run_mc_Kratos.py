import sys
sys.dont_write_bytecode = True
import os
import json

import xmc

if __name__ == "__main__":

    if(len(sys.argv)==2):
        parametersPath = str(sys.argv[1]) # set path to the parameters
    else:
        parametersPath = "parameters/parameters_xmc_test_mc_Kratos_asynchronous_poisson_2d.json"

    # read parameters
    with open(parametersPath,'r') as parameter_file:
            parameters = json.load(parameter_file)

    # add path of the problem folder to python path
    problem_id = parameters["solverWrapperInputDictionary"]["problemId"]
    sys.path.append(os.path.join("..","xmc","classDefs_solverWrapper","problemDefs_KratosMultiphysics",problem_id))

    # SampleGenerator
    samplerInputDictionary = parameters["samplerInputDictionary"]
    samplerInputDictionary['randomGeneratorInputDictionary'] = parameters["randomGeneratorInputDictionary"]
    samplerInputDictionary['solverWrapperInputDictionary'] = parameters["solverWrapperInputDictionary"]

    # MonteCarloIndex
    monteCarloIndexInputDictionary = parameters["monteCarloIndexInputDictionary"]
    monteCarloIndexInputDictionary["samplerInputDictionary"] = samplerInputDictionary

    # Moment Estimators
    qoiEstimatorInputDictionary = parameters["qoiEstimatorInputDictionary"]
    combinedEstimatorInputDictionary = parameters["combinedEstimatorInputDictionary"]
    costEstimatorInputDictionary = parameters["costEstimatorInputDictionary"]
    # qoi estimators
    monteCarloIndexInputDictionary["qoiEstimator"] = [monteCarloIndexInputDictionary["qoiEstimator"][0] for _ in range (0,parameters["solverWrapperInputDictionary"]["numberQoI"])]
    monteCarloIndexInputDictionary["qoiEstimatorInputDictionary"] = [qoiEstimatorInputDictionary]*parameters["solverWrapperInputDictionary"]["numberQoI"]
    # combined estimators
    monteCarloIndexInputDictionary["combinedEstimator"] = [monteCarloIndexInputDictionary["combinedEstimator"][0] for _ in range (0,parameters["solverWrapperInputDictionary"]["numberCombinedQoi"])]
    monteCarloIndexInputDictionary["combinedEstimatorInputDictionary"] = [combinedEstimatorInputDictionary]*parameters["solverWrapperInputDictionary"]["numberCombinedQoi"]
    # cost estimator
    monteCarloIndexInputDictionary["costEstimatorInputDictionary"] = costEstimatorInputDictionary

    # MonoCriterion
    criteriaArray = []
    criteriaInputs = []
    for monoCriterion in (parameters["monoCriteriaInpuctDictionary"]):
        criteriaArray.append(xmc.monoCriterion.MonoCriterion(\
            parameters["monoCriteriaInpuctDictionary"][monoCriterion]["criteria"],\
            parameters["monoCriteriaInpuctDictionary"][monoCriterion]["tolerance"]))
        criteriaInputs.append([parameters["monoCriteriaInpuctDictionary"][monoCriterion]["input"]])

    # MultiCriterion
    multiCriterionInputDictionary=parameters["multiCriterionInputDictionary"]
    multiCriterionInputDictionary["criteria"] = criteriaArray
    multiCriterionInputDictionary["inputsForCriterion"] = criteriaInputs
    criterion = xmc.multiCriterion.MultiCriterion(**multiCriterionInputDictionary)

    # ErrorEstimator
    statErrorEstimator = xmc.errorEstimator.ErrorEstimator(**parameters["errorEstimatorInputDictionary"])

    # HierarchyOptimiser
    hierarchyCostOptimiser = xmc.hierarchyOptimiser.HierarchyOptimiser(**parameters["hierarchyOptimiserInputDictionary"])

    # EstimationAssembler
    if "expectationAssembler" in parameters["estimationAssemblerInputDictionary"].keys():
        expectationAssembler = xmc.estimationAssembler.EstimationAssembler(**parameters["estimationAssemblerInputDictionary"]["expectationAssembler"])
    if "varianceAssembler" in parameters["estimationAssemblerInputDictionary"].keys():
        varianceAssembler = xmc.estimationAssembler.EstimationAssembler(**parameters["estimationAssemblerInputDictionary"]["varianceAssembler"])

    # MonteCarloSampler
    monteCarloSamplerInputDictionary = parameters["monteCarloSamplerInputDictionary"]
    monteCarloSamplerInputDictionary["indexConstructorDictionary"] = monteCarloIndexInputDictionary
    monteCarloSamplerInputDictionary["assemblers"] =  [expectationAssembler,varianceAssembler]
    monteCarloSamplerInputDictionary["errorEstimators"] = [statErrorEstimator]
    mcSampler = xmc.monteCarloSampler.MonteCarloSampler(**monteCarloSamplerInputDictionary)

    # XMCAlgorithm
    XMCAlgorithmInputDictionary = parameters["XMCAlgorithmInputDictionary"]
    XMCAlgorithmInputDictionary["monteCarloSampler"] = mcSampler
    XMCAlgorithmInputDictionary["hierarchyOptimiser"] = hierarchyCostOptimiser
    XMCAlgorithmInputDictionary["stoppingCriterion"] = criterion

    algo = xmc.XMCAlgorithm(**XMCAlgorithmInputDictionary)

    if (parameters["solverWrapperInputDictionary"]["asynchronous"] is True):
        algo.runAsynchronousXMC()
    else:
        algo.runXMC()

    ########################################################################################################################################################################################################
    ########################################################################################################################################################################################################
    ########################################################################################################################################################################################################

    # from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
    # import xmc.methodDefs_momentEstimator.computeCentralMoments as mdccm

    # # writing to file a dictionary
    # qoi_dict = {}
    # for qoi_counter in range (parameters["solverWrapperInputDictionary"]["numberQoI"]):
    #     qoi_dict[qoi_counter] = {index: {} for index in range (len(algo.monteCarloSampler.indices))}
    #     for index in range (len(algo.monteCarloSampler.indices)):
    #         algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter] = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter])
    #         sample_counter = algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._sampleCounter
    #         S1 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[0][0])
    #         S2 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[1][0])
    #         h1 = get_value_from_remote(mdccm.computeCentralMomentsOrderOneDimensionZero(S1,sample_counter))
    #         h2 = get_value_from_remote(mdccm.computeCentralMomentsOrderTwoDimensionZero(S1,S2,sample_counter))
    #         qoi_dict[qoi_counter][index] = {"qoi_id":qoi_counter, "index": index, "instances": sample_counter, "S1": S1, "S2": S2, "h1": h1, "h2": h2}
    # for qoi_counter in range (parameters["solverWrapperInputDictionary"]["numberQoI"],parameters["solverWrapperInputDictionary"]["numberQoI"]+parameters["solverWrapperInputDictionary"]["numberCombinedQoi"]):
    #         qoi_dict[qoi_counter] = {index: {} for index in range (len(algo.monteCarloSampler.indices))}
    #         for index in range (len(algo.monteCarloSampler.indices)):
    #             algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter] = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter])
    #             sample_counter = algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._sampleCounter
    #             S1 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[0][0])
    #             S2 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[1][0])
    #             h1 = get_value_from_remote(mdccm.computeCentralMomentsOrderOneDimensionZero(S1,sample_counter))
    #             h2 = get_value_from_remote(mdccm.computeCentralMomentsOrderTwoDimensionZeroBiased(S1,S2,sample_counter))
    #             qoi_dict[qoi_counter][index] = {"qoi_id":qoi_counter, "index": index, "instances": sample_counter, "S1": S1, "S2": S2, "h1": h1, "h2": h2}

    # with open('MC_asynchronous_power_sums.json', 'w') as f:
    #     json.dump(qoi_dict, f, indent=2)