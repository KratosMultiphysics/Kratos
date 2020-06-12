# Import python class test
import unittest

# Import python libraries
import json
import sys
import os

# Import xmc classes
import xmc

# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import get_value_from_remote   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import get_value_from_remote  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import get_value_from_remote      # to execute with python3


class TestXMCAlgorithm(unittest.TestCase):

    def test_mc_asynchronous_Kratos(self):

        # read parameters
        parametersList = ["parameters/parameters_xmc_test_mc_Kratos_asynchronous_poisson_2d.json", \
            "parameters/parameters_xmc_test_mc_Kratos_asynchronous_poisson_2d_with_combined_power_sums.json", \
            "parameters/parameters_xmc_test_mc_Kratos_asynchronous_poisson_2d_with_10_combined_power_sums.json", \
            "parameters/parameters_xmc_test_mc_Kratos_poisson_2d.json", \
            "parameters/parameters_xmc_test_mc_Kratos_poisson_2d_with_combined_power_sums.json"]

        for parametersPath in parametersList:
            with open(parametersPath,'r') as parameter_file:
                parameters = json.load(parameter_file)
            # add path of the problem folder to python path
            problem_id = parameters["solverWrapperInputDictionary"]["problemId"]
            sys.path.append(os.path.join("poisson_square_2d_xmc"))
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
            for monoCriterion in (parameters["monoCriteriaInpuctDict"]):
                criteriaArray.append(xmc.monoCriterion.MonoCriterion(\
                    parameters["monoCriteriaInpuctDict"][monoCriterion]["criteria"],\
                    parameters["monoCriteriaInpuctDict"][monoCriterion]["tolerance"]))
                criteriaInputs.append([parameters["monoCriteriaInpuctDict"][monoCriterion]["input"]])
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

            # test
            estimations = get_value_from_remote(algo.estimation())
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
            problem_id = parameters["solverWrapperInputDictionary"]["problemId"]
            sys.path.append(os.path.join("poisson_square_2d_xmc"))
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
            for monoCriterion in (parameters["monoCriteriaInpuctDict"]):
                criteriaArray.append(xmc.monoCriterion.MonoCriterion(\
                    parameters["monoCriteriaInpuctDict"][monoCriterion]["criteria"],\
                    parameters["monoCriteriaInpuctDict"][monoCriterion]["tolerance"]))
                criteriaInputs.append([parameters["monoCriteriaInpuctDict"][monoCriterion]["input"]])
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

            # test
            estimations = get_value_from_remote(algo.estimation())
            estimated_mean = 1.47
            self.assertAlmostEqual(estimations[0],estimated_mean,delta=1.0)
            self.assertEqual(algo.monteCarloSampler.indices[0].costEstimator._sampleCounter,15) # level 0
            self.assertEqual(algo.monteCarloSampler.indices[1].costEstimator._sampleCounter,15) # level 1
            self.assertEqual(algo.monteCarloSampler.indices[2].costEstimator._sampleCounter,15) # level 2

if __name__ == '__main__':
    unittest.main()
