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
    problem_id = parameters["solverWrapperInputDict"]["problemId"]
    sys.path.append(os.path.join("..","xmc","classDefs_solverWrapper","problemDefs_KratosMultiphysics",problem_id))

    # RandomGeneratorWrapper
    randomGeneratorInputDict = parameters["randomGeneratorInputDict"]

    solverWrapperInputDict = parameters["solverWrapperInputDict"]

    # SampleGenerator
    samplerInputDict = parameters["samplerInputDict"]
    samplerInputDict['randomGeneratorInputDict'] = randomGeneratorInputDict
    samplerInputDict['solverWrapperInputDict'] = solverWrapperInputDict

    # Moment Estimators
    qoiEstimatorInputDict = parameters["qoiEstimatorInputDict"]
    combinedEstimatorInputDict = parameters["combinedEstimatorInputDict"]
    costEstimatorInputDict = parameters["costEstimatorInputDict"]

    # MonteCarloIndex Constructor
    monteCarloIndexInputDict = parameters["monteCarloIndexInputDict"]
    monteCarloIndexInputDict["samplerInputDict"] = samplerInputDict
    # qoi estimators
    monteCarloIndexInputDict["qoiEstimator"] = [monteCarloIndexInputDict["qoiEstimator"][0] for _ in range (0,parameters["solverWrapperInputDict"]["numberQoI"])]
    monteCarloIndexInputDict["qoiEstimatorInputDict"] = [qoiEstimatorInputDict]*parameters["solverWrapperInputDict"]["numberQoI"]
    # combined estimators
    monteCarloIndexInputDict["combinedEstimator"] = [monteCarloIndexInputDict["combinedEstimator"][0] for _ in range (0,parameters["solverWrapperInputDict"]["numberCombinedQoi"])]
    monteCarloIndexInputDict["combinedEstimatorInputDict"] = [combinedEstimatorInputDict]*parameters["solverWrapperInputDict"]["numberCombinedQoi"]
    # cost estimator
    monteCarloIndexInputDict["costEstimatorInputDict"] = costEstimatorInputDict

    ################################# RUN TIME GENERATED ENTITIES END HERE #######################

    # MonoCriterion
    criteriaArray = []
    criteriaInputs = []
    for monoCriterion in (parameters["monoCriteriaInpuctDict"]):
        criteriaArray.append(xmc.monoCriterion.MonoCriterion(\
            parameters["monoCriteriaInpuctDict"][monoCriterion]["criteria"],\
            parameters["monoCriteriaInpuctDict"][monoCriterion]["tolerance"]))
        criteriaInputs.append([parameters["monoCriteriaInpuctDict"][monoCriterion]["input"]])

    # MultiCriterion
    criterion = xmc.multiCriterion.MultiCriterion(criteria=criteriaArray,
                                                  inputsForCriterion=criteriaInputs,
                                                  interpreter='xmc.methodDefs_multiCriterion.interpreter.interpretAsConvergenceAndIterationBounds',
                                                  flag='xmc.methodDefs_multiCriterion.flag.plainFlag')


    # ErrorEstimator
    # TODO we should define xmc.methodDefs_errorEstimator.errorEstimation.Variance+BiasError_Task
    MSEErrorEstimator = xmc.errorEstimator.ErrorEstimator(
        error='xmc.methodDefs_errorEstimator.errorEstimation.errorEstimationMSE_Task',
        parameters=[0.95])

    # HierarchyOptimiser
    hierarchyOptimiserInputDict = parameters["hierarchyOptimiserInputDict"]
    hierarchyCostOptimiser = xmc.hierarchyOptimiser.HierarchyOptimiser(**hierarchyOptimiserInputDict)


    # EstimationAssembler
    expectationAssembler = xmc.estimationAssembler.EstimationAssembler(
        assembleEstimation='xmc.methodDefs_estimationAssembler.assembleEstimation.assembleValue_Task')
    biasAssembler = xmc.estimationAssembler.EstimationAssembler(
        assembleEstimation='xmc.methodDefs_estimationAssembler.assembleEstimation.assembleBias_Task')
    varianceAssembler = xmc.estimationAssembler.EstimationAssembler(
        assembleEstimation='xmc.methodDefs_estimationAssembler.assembleEstimation.assembleStatisticalError_Task')

    # MonteCarloSampler
    monteCarloSamplerInputDict = parameters["monteCarloSamplerInputDict"]
    monteCarloSamplerInputDict["indexConstructorDictionary"] = monteCarloIndexInputDict
    monteCarloSamplerInputDict["assemblers"] =  [expectationAssembler,biasAssembler,varianceAssembler]
    monteCarloSamplerInputDict["errorEstimators"] = [MSEErrorEstimator]
    mcSampler = xmc.monteCarloSampler.MonteCarloSampler(**monteCarloSamplerInputDict)

    # XMCAlgorithm
    XMCAlgorithmInputDict = parameters["XMCAlgorithmInputDict"]
    XMCAlgorithmInputDict["monteCarloSampler"] = mcSampler
    XMCAlgorithmInputDict["hierarchyOptimiser"] = hierarchyCostOptimiser
    XMCAlgorithmInputDict["stoppingCriterion"] = criterion

    algo = xmc.XMCAlgorithm(**XMCAlgorithmInputDict)

    if (parameters["solverWrapperInputDict"]["asynchronous"] is True):
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
    # for qoi_counter in range (parameters["solverWrapperInputDict"]["numberQoI"]):
    #     qoi_dict[qoi_counter] = {index: {} for index in range (len(algo.monteCarloSampler.indices))}
    #     for index in range (len(algo.monteCarloSampler.indices)):
    #         sample_counter = algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._sampleCounter
    #         S10 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[0][0])
    #         S01 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[0][1])
    #         S20 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[1][0])
    #         S11 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[1][1])
    #         S02 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[1][2])
    #         h1 = get_value_from_remote(mdccm.computeCentralMomentsOrderOneDimensionZero_Task(S01,sample_counter))
    #         h2 = get_value_from_remote(mdccm.computeCentralMomentsOrderTwoDimensionZero_Task(S01,S02,sample_counter))
    #         qoi_dict[qoi_counter][index] = {"qoi_id":qoi_counter, "index": index, "instances": sample_counter, "S10": S10, "S01": S01, "S20": S20, "S11": S11, "S02": S02, "h1": h1, "h2": h2}
    # for qoi_counter in range (parameters["solverWrapperInputDict"]["numberQoI"],parameters["solverWrapperInputDict"]["numberQoI"]+parameters["solverWrapperInputDict"]["numberTimePowerSums"]):
    #     qoi_dict[qoi_counter] = {index: {} for index in range (len(algo.monteCarloSampler.indices))}
    #     for index in range (len(algo.monteCarloSampler.indices)):
    #         if (type(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._sampleCounter) is list):
    #             sample_counter = 0
    #             for i in range (len(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._sampleCounter)):
    #                 sample_counter = sample_counter + get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._sampleCounter[i])
    #         else:
    #             sample_counter = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._sampleCounter)
    #         S1upper = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[0][0])
    #         S2upper = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[1][0])
    #         S1lower = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[0][1])
    #         S2lower = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[1][1])
    #         h1upper = get_value_from_remote(mdccm.computeCentralMomentsOrderOneDimensionZero_Task(S1upper,sample_counter))
    #         h2upper = get_value_from_remote(mdccm.computeCentralMomentsOrderTwoDimensionZeroBiased_Task(S1upper,S2upper,sample_counter))
    #         h1lower = get_value_from_remote(mdccm.computeCentralMomentsOrderOneDimensionZero_Task(S1lower,sample_counter))
    #         h2lower = get_value_from_remote(mdccm.computeCentralMomentsOrderTwoDimensionZeroBiased_Task(S1lower,S2lower,sample_counter))
    #         qoi_dict[qoi_counter][index] = {"qoi_id":qoi_counter, "index": index, "instances": sample_counter, "S1_upper": S1upper, "S1_lower": S1lower, "S2_upper": S2upper, "S2_lower": S2lower, "h1_upper": h1upper, "h1_lower": h1lower, "h2_upper": h2upper, "h2_lower": h2lower, "h1": (h1upper-h1lower), "h2": (h2upper-h2lower)}
    # with open('MLMC_asynchronous_power_sums.json', 'w') as f:
    #     json.dump(qoi_dict, f, indent=2)