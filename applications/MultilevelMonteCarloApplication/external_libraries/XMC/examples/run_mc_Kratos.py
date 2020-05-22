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

    #################### RUN TIME GENERATED ENTITIES END HERE ####################

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
    statErrorEstimator = xmc.errorEstimator.ErrorEstimator(
        error='xmc.methodDefs_errorEstimator.errorEstimation.errorEstimationStatError_Task',
        parameters=[0.95])

    # HierarchyOptimiser
    hierarchyOptimiserInputDict = parameters["hierarchyOptimiserInputDict"]
    hierarchyCostOptimiser = xmc.hierarchyOptimiser.HierarchyOptimiser(**hierarchyOptimiserInputDict)

    # EstimationAssembler
    expectationAssembler = xmc.estimationAssembler.EstimationAssembler(
        assembleEstimation='xmc.methodDefs_estimationAssembler.assembleEstimation.assembleValue_Task')
    varianceAssembler = xmc.estimationAssembler.EstimationAssembler(
        assembleEstimation='xmc.methodDefs_estimationAssembler.assembleEstimation.assembleStatisticalError_Task')

    # MonteCarloSampler
    monteCarloSamplerInputDict = parameters["monteCarloSamplerInputDict"]
    monteCarloSamplerInputDict["indexConstructorDictionary"] = monteCarloIndexInputDict
    monteCarloSamplerInputDict["assemblers"] =  [expectationAssembler,varianceAssembler]
    monteCarloSamplerInputDict["errorEstimators"] = [statErrorEstimator]
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
    #         S1 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[0][0])
    #         S2 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[1][0])
    #         h1 = get_value_from_remote(mdccm.computeCentralMomentsOrderOneDimensionZero_Task(S1,sample_counter))
    #         h2 = get_value_from_remote(mdccm.computeCentralMomentsOrderTwoDimensionZero_Task(S1,S2,sample_counter))
    #         qoi_dict[qoi_counter][index] = {"qoi_id":qoi_counter, "index": index, "instances": sample_counter, "S1": S1, "S2": S2, "h1": h1, "h2": h2}
    # for qoi_counter in range (parameters["solverWrapperInputDict"]["numberQoI"],parameters["solverWrapperInputDict"]["numberQoI"]+parameters["solverWrapperInputDict"]["numberCombinedQoi"]):
    #         qoi_dict[qoi_counter] = {index: {} for index in range (len(algo.monteCarloSampler.indices))}
    #         for index in range (len(algo.monteCarloSampler.indices)):
    #             if (type(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._sampleCounter) is list):
    #                 sample_counter = 0
    #                 for i in range (len(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._sampleCounter)):
    #                     sample_counter = sample_counter + get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._sampleCounter[i])
    #             else:
    #                 sample_counter = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter]._sampleCounter)
    #             S1 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[0][0])
    #             S2 = get_value_from_remote(algo.monteCarloSampler.indices[index].qoiEstimator[qoi_counter].powerSums[1][0])
    #             h1 = get_value_from_remote(mdccm.computeCentralMomentsOrderOneDimensionZero_Task(S1,sample_counter))
    #             h2 = get_value_from_remote(mdccm.computeCentralMomentsOrderTwoDimensionZeroBiased_Task(S1,S2,sample_counter))
    #             qoi_dict[qoi_counter][index] = {"qoi_id":qoi_counter, "index": index, "instances": sample_counter, "S1": S1, "S2": S2, "h1": h1, "h2": h2}
    # with open('MC_asynchronous_power_sums.json', 'w') as f:
    #     json.dump(qoi_dict, f, indent=2)