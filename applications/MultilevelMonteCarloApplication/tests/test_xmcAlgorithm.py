# Import python class test
import unittest

# Import python libraries
import json
import sys
import os

# Import xmc classes
import xmc


class TestXMCAlgorithm(unittest.TestCase):

    def test_mc_asynchronous_Kratos(self):

        # read parameters
        parametersList = ["parameters/parameters_xmc_test_mc_Kratos_asynchronous_poisson_2d.json", \
            "parameters/parameters_xmc_test_mc_Kratos_asynchronous_poisson_2d_with_combined_power_sums.json", \
            "parameters/parameters_xmc_test_mc_Kratos_poisson_2d.json", \
            "parameters/parameters_xmc_test_mc_Kratos_poisson_2d_with_combined_power_sums.json"]

        for parametersPath in parametersList:
            with open(parametersPath,'r') as parameter_file:
                parameters = json.load(parameter_file)
            # add path of the problem folder to python path
            problem_id = parameters["solverWrapperInputDict"]["problemId"]
            sys.path.append(os.path.join("poisson_square_2d_xmc"))
            # RandomGeneratorWrapper
            randomGeneratorInputDict = parameters["randomGeneratorInputDict"]
            # SolverWrapper
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
            # run
            if (parameters["solverWrapperInputDict"]["asynchronous"] is True):
                algo.runAsynchronousXMC()
            else:
                algo.runXMC()

            # test
            estimations = algo.estimation()
            estimated_mean = 1.5
            self.assertAlmostEqual(estimations[0],estimated_mean,delta=0.1)
            self.assertEqual(algo.monteCarloSampler.indices[0].costEstimator._sampleCounter,15)


    def test_mlmc_asynchronous_Kratos(self):

        # read parameters
        parametersList = ["parameters/parameters_xmc_test_mlmc_Kratos_asynchronous_poisson_2d.json", \
        "parameters/parameters_xmc_test_mlmc_Kratos_asynchronous_poisson_2d_with_combined_power_sums.json", \
        "parameters/parameters_xmc_test_mlmc_Kratos_poisson_2d.json", \
        "parameters/parameters_xmc_test_mlmc_Kratos_poisson_2d_with_combined_power_sums.json"]
        for parametersPath in parametersList:
            with open(parametersPath,'r') as parameter_file:
                    parameters = json.load(parameter_file)
            # add path of the problem folder to python path
            problem_id = parameters["solverWrapperInputDict"]["problemId"]
            sys.path.append(os.path.join("poisson_square_2d_xmc"))
            # RandomGeneratorWrapper
            randomGeneratorInputDict = parameters["randomGeneratorInputDict"]
            # SolverWrapper
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
            # run
            if (parameters["solverWrapperInputDict"]["asynchronous"] is True):
                algo.runAsynchronousXMC()
            else:
                algo.runXMC()

            # test
            estimations = algo.estimation()
            estimated_mean = 1.47
            self.assertAlmostEqual(estimations[0],estimated_mean,delta=1.0)
            self.assertEqual(algo.monteCarloSampler.indices[0].costEstimator._sampleCounter,15) # level 0
            self.assertEqual(algo.monteCarloSampler.indices[1].costEstimator._sampleCounter,15) # level 1
            self.assertEqual(algo.monteCarloSampler.indices[2].costEstimator._sampleCounter,15) # level 2

if __name__ == '__main__':
    unittest.main()
