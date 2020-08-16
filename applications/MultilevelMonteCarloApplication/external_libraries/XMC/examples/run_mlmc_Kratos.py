import sys
sys.dont_write_bytecode = True
import os
import json

import xmc

if __name__ == "__main__":

    if(len(sys.argv)==2):
        parametersPath = str(sys.argv[1]) # set path to the parameters
    else:
        parametersPath = "parameters/parameters_xmc_test_mlmc_Kratos_asynchronous_poisson_2d.json"

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

    # MonteCarloIndex Constructor
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
    MSEErrorEstimator = xmc.errorEstimator.ErrorEstimator(**parameters["errorEstimatorInputDictionary"])

    # HierarchyOptimiser
    hierarchyCostOptimiser = xmc.hierarchyOptimiser.HierarchyOptimiser(**parameters["hierarchyOptimiserInputDictionary"])

    # EstimationAssembler
    if "expectationAssembler" in parameters["estimationAssemblerInputDictionary"].keys():
        expectationAssembler = xmc.estimationAssembler.EstimationAssembler(**parameters["estimationAssemblerInputDictionary"]["expectationAssembler"])
    if "discretizationErrorAssembler" in parameters["estimationAssemblerInputDictionary"].keys():
        discretizationErrorAssembler = xmc.estimationAssembler.EstimationAssembler(**parameters["estimationAssemblerInputDictionary"]["discretizationErrorAssembler"])
    if "varianceAssembler" in parameters["estimationAssemblerInputDictionary"].keys():
        varianceAssembler = xmc.estimationAssembler.EstimationAssembler(**parameters["estimationAssemblerInputDictionary"]["varianceAssembler"])

    # MonteCarloSampler
    monteCarloSamplerInputDictionary = parameters["monteCarloSamplerInputDictionary"]
    monteCarloSamplerInputDictionary["indexConstructorDictionary"] = monteCarloIndexInputDictionary
    monteCarloSamplerInputDictionary["assemblers"] =  [expectationAssembler,discretizationErrorAssembler,varianceAssembler]
    monteCarloSamplerInputDictionary["errorEstimators"] = [MSEErrorEstimator]
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
    #         S10 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[0][0])
    #         S01 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[0][1])
    #         S20 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[1][0])
    #         S11 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[1][1])
    #         S02 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[1][2])
    #         h1 = get_value_from_remote(mdccm.computeCentralMomentsOrderOneDimensionZero(S01,sample_counter))
    #         h2 = get_value_from_remote(mdccm.computeCentralMomentsOrderTwoDimensionZero(S01,S02,sample_counter))
    #         qoi_dict[qoi_counter][index] = {"qoi_id":qoi_counter, "index": index, "instances": sample_counter, "S10": S10, "S01": S01, "S20": S20, "S11": S11, "S02": S02, "h1": h1, "h2": h2}
    # for qoi_counter in range (parameters["solverWrapperInputDictionary"]["numberQoI"],parameters["solverWrapperInputDictionary"]["numberQoI"]+parameters["solverWrapperInputDictionary"]["numberCombinedQoi"]):
    #     qoi_dict[qoi_counter] = {index: {} for index in range (len(algo.monteCarloSampler.indices))}
    #     for index in range (len(algo.monteCarloSampler.indices)):
    #         algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter] = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter])
    #         sample_counter = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._sampleCounter)
    #         S1upper = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[0][0])
    #         S2upper = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[1][0])
    #         S1lower = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[0][1])
    #         S2lower = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[1][1])
    #         h1upper = get_value_from_remote(mdccm.computeCentralMomentsOrderOneDimensionZero(S1upper,sample_counter))
    #         h2upper = get_value_from_remote(mdccm.computeCentralMomentsOrderTwoDimensionZeroBiased(S1upper,S2upper,sample_counter))
    #         h1lower = get_value_from_remote(mdccm.computeCentralMomentsOrderOneDimensionZero(S1lower,sample_counter))
    #         h2lower = get_value_from_remote(mdccm.computeCentralMomentsOrderTwoDimensionZeroBiased(S1lower,S2lower,sample_counter))
    #         qoi_dict[qoi_counter][index] = {"qoi_id":qoi_counter, "index": index, "instances": sample_counter, "S1_upper": S1upper, "S1_lower": S1lower, "S2_upper": S2upper, "S2_lower": S2lower, "h1_upper": h1upper, "h1_lower": h1lower, "h2_upper": h2upper, "h2_lower": h2lower, "h1": (h1upper-h1lower), "h2": (h2upper-h2lower)}
    # with open('MLMC_asynchronous_power_sums.json', 'w') as f:
    #     json.dump(qoi_dict, f, indent=2)