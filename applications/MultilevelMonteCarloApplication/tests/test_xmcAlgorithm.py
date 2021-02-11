# Import Python libraries
import unittest
import json
import sys
import os

# Import XMC, distributed environment
import xmc
from xmc.distributedEnvironmentFramework import get_value_from_remote

def kratosFound():
    try:
        from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
        return CheckIfApplicationsAvailable(
            "ConvectionDiffusionApplication",
            "LinearSolversApplication",
            "MeshingApplication",
            "MultilevelMonteCarloApplication",
        )
    except ImportError:
        return False

class TestXMCAlgorithm(unittest.TestCase):
    def test_mc_Kratos(self):
        if not kratosFound():
            self.skipTest("Missing dependency: KratosMultiphysics")

        # read parameters
        parametersList = [
            "parameters/parameters_xmc_test_mc_Kratos_asynchronous_poisson_2d.json",
            "parameters/parameters_xmc_test_mc_Kratos_asynchronous_poisson_2d_with_combined_power_sums.json",
            "parameters/parameters_xmc_test_mc_Kratos_asynchronous_poisson_2d_with_10_combined_power_sums.json",
            "parameters/parameters_xmc_test_mc_Kratos_poisson_2d.json",
            "parameters/parameters_xmc_test_mc_Kratos_poisson_2d_with_combined_power_sums.json",
            "parameters/poisson_multi-moment_mc.json",
        ]

        for parametersPath in parametersList:
            with open(parametersPath, "r") as parameter_file:
                parameters = json.load(parameter_file)
            # add path of the problem folder to python path
            problem_id = parameters["solverWrapperInputDictionary"]["problemId"]
            sys.path.append(os.path.join("poisson_square_2d"))
            # SampleGenerator
            samplerInputDictionary = parameters["samplerInputDictionary"]
            samplerInputDictionary["randomGeneratorInputDictionary"] = parameters[
                "randomGeneratorInputDictionary"
            ]
            samplerInputDictionary["solverWrapperInputDictionary"] = parameters[
                "solverWrapperInputDictionary"
            ]
            # MonteCarloIndex
            monteCarloIndexInputDictionary = parameters["monteCarloIndexInputDictionary"]
            monteCarloIndexInputDictionary["samplerInputDictionary"] = samplerInputDictionary
            # Moment Estimators
            if "qoiEstimatorInputDictionary" not in monteCarloIndexInputDictionary.keys():
                numberMomentEstimator = parameters["solverWrapperInputDictionary"]["numberMomentEstimator"]
                monteCarloIndexInputDictionary["qoiEstimator"] = [
                    monteCarloIndexInputDictionary["qoiEstimator"][0] for _ in range(numberMomentEstimator)
                ]
                monteCarloIndexInputDictionary["qoiEstimatorInputDictionary"] = [
                    parameters["qoiEstimatorInputDictionary"]
                ] * numberMomentEstimator
            # combined estimators
            if "combinedEstimator" in monteCarloIndexInputDictionary.keys():
                numberCombinedQoI = parameters["solverWrapperInputDictionary"][
                    "numberCombinedMomentEstimator"
                ]
                monteCarloIndexInputDictionary["combinedEstimator"] = [
                    monteCarloIndexInputDictionary["combinedEstimator"][0]
                    for _ in range(numberCombinedQoI)
                ]
                monteCarloIndexInputDictionary["combinedEstimatorInputDictionary"] = [
                    parameters["combinedEstimatorInputDictionary"]
                ] * numberCombinedQoI
            # cost estimator
            monteCarloIndexInputDictionary["costEstimatorInputDictionary"] = parameters[
                "costEstimatorInputDictionary"
            ]
            # MonoCriterion
            criteriaArray = []
            criteriaInputs = []
            for monoCriterion in parameters["monoCriteriaInpuctDictionary"]:
                criteriaArray.append(
                    xmc.monoCriterion.MonoCriterion(
                        parameters["monoCriteriaInpuctDictionary"][monoCriterion]["criteria"],
                        parameters["monoCriteriaInpuctDictionary"][monoCriterion]["tolerance"],
                    )
                )
                criteriaInputs.append(
                    [parameters["monoCriteriaInpuctDictionary"][monoCriterion]["input"]]
                )
            # MultiCriterion
            multiCriterionInputDictionary = parameters["multiCriterionInputDictionary"]
            multiCriterionInputDictionary["criteria"] = criteriaArray
            multiCriterionInputDictionary["inputsForCriterion"] = criteriaInputs
            criterion = xmc.multiCriterion.MultiCriterion(**multiCriterionInputDictionary)
            # ErrorEstimator
            statErrorEstimator = xmc.errorEstimator.ErrorEstimator(
                **parameters["errorEstimatorInputDictionary"]
            )
            # HierarchyOptimiser
            hierarchyCostOptimiser = xmc.hierarchyOptimiser.HierarchyOptimiser(
                **parameters["hierarchyOptimiserInputDictionary"]
            )
            # EstimationAssembler
            if (
                "expectationAssembler"
                in parameters["estimationAssemblerInputDictionary"].keys()
            ):
                expectationAssembler = xmc.estimationAssembler.EstimationAssembler(
                    **parameters["estimationAssemblerInputDictionary"]["expectationAssembler"]
                )
            if "varianceAssembler" in parameters["estimationAssemblerInputDictionary"].keys():
                varianceAssembler = xmc.estimationAssembler.EstimationAssembler(
                    **parameters["estimationAssemblerInputDictionary"]["varianceAssembler"]
                )
            # MonteCarloSampler
            monteCarloSamplerInputDictionary = parameters["monteCarloSamplerInputDictionary"]
            monteCarloSamplerInputDictionary[
                "indexConstructorDictionary"
            ] = monteCarloIndexInputDictionary
            monteCarloSamplerInputDictionary["assemblers"] = [
                expectationAssembler,
                varianceAssembler,
            ]
            monteCarloSamplerInputDictionary["errorEstimators"] = [statErrorEstimator]
            mcSampler = xmc.monteCarloSampler.MonteCarloSampler(
                **monteCarloSamplerInputDictionary
            )
            # XMCAlgorithm
            XMCAlgorithmInputDictionary = parameters["XMCAlgorithmInputDictionary"]
            XMCAlgorithmInputDictionary["monteCarloSampler"] = mcSampler
            XMCAlgorithmInputDictionary["hierarchyOptimiser"] = hierarchyCostOptimiser
            XMCAlgorithmInputDictionary["stoppingCriterion"] = criterion
            algo = xmc.XMCAlgorithm(**XMCAlgorithmInputDictionary)

            if parameters["solverWrapperInputDictionary"]["asynchronous"] is True:
                algo.runAsynchronousXMC()
            else:
                algo.runXMC()

            # test
            estimations = get_value_from_remote(algo.estimation())
            estimated_mean = 1.5
            self.assertAlmostEqual(estimations[0], estimated_mean, delta=0.1)
            self.assertEqual(algo.hierarchy()[0][1], 15)

    def test_mlmc_Kratos(self):
        if not kratosFound():
            self.skipTest("Missing dependency: KratosMultiphysics")

        # read parameters
        parametersList = [
            "parameters/parameters_xmc_test_mlmc_Kratos_asynchronous_poisson_2d.json",
            "parameters/parameters_xmc_test_mlmc_Kratos_asynchronous_poisson_2d_with_combined_power_sums.json",
            "parameters/parameters_xmc_test_mlmc_Kratos_poisson_2d.json",
            "parameters/parameters_xmc_test_mlmc_Kratos_poisson_2d_with_combined_power_sums.json",
            "parameters/poisson_multi-moment_mlmc.json",
            "parameters/parameters_xmc_test_mlmc_Kratos_asynchronous_poisson_2d_with_combined_power_sums_multi.json",
            "parameters/parameters_xmc_test_mlmc_Kratos_asynchronous_poisson_2d_DAR.json"
        ]
        for parametersPath in parametersList:
            with open(parametersPath, "r") as parameter_file:
                parameters = json.load(parameter_file)
            # add path of the problem folder to python path
            problem_id = parameters["solverWrapperInputDictionary"]["problemId"]
            sys.path.append(os.path.join("poisson_square_2d"))
            # SampleGenerator
            samplerInputDictionary = parameters["samplerInputDictionary"]
            samplerInputDictionary["randomGeneratorInputDictionary"] = parameters[
                "randomGeneratorInputDictionary"
            ]
            samplerInputDictionary["solverWrapperInputDictionary"] = parameters[
                "solverWrapperInputDictionary"
            ]
            # MonteCarloIndex Constructor
            monteCarloIndexInputDictionary = parameters["monteCarloIndexInputDictionary"]
            monteCarloIndexInputDictionary["samplerInputDictionary"] = samplerInputDictionary
            # Moment Estimators
            if "qoiEstimatorInputDictionary" not in monteCarloIndexInputDictionary.keys():
                numberMomentEstimator = parameters["solverWrapperInputDictionary"]["numberMomentEstimator"]
                monteCarloIndexInputDictionary["qoiEstimator"] = [
                    monteCarloIndexInputDictionary["qoiEstimator"][0] for _ in range(numberMomentEstimator)
                ]
                monteCarloIndexInputDictionary["qoiEstimatorInputDictionary"] = [
                    parameters["qoiEstimatorInputDictionary"]
                ] * numberMomentEstimator
            # combined estimators
            if "combinedEstimator" in monteCarloIndexInputDictionary.keys():
                numberCombinedQoI = parameters["solverWrapperInputDictionary"][
                    "numberCombinedMomentEstimator"
                ]
                monteCarloIndexInputDictionary["combinedEstimator"] = [
                    monteCarloIndexInputDictionary["combinedEstimator"][0]
                    for _ in range(numberCombinedQoI)
                ]
                monteCarloIndexInputDictionary["combinedEstimatorInputDictionary"] = [
                    parameters["combinedEstimatorInputDictionary"]
                ] * numberCombinedQoI
            # cost estimator
            monteCarloIndexInputDictionary["costEstimatorInputDictionary"] = parameters[
                "costEstimatorInputDictionary"
            ]
            # MonoCriterion
            criteriaArray = []
            criteriaInputs = []
            for monoCriterion in parameters["monoCriteriaInpuctDictionary"]:
                criteriaArray.append(
                    xmc.monoCriterion.MonoCriterion(
                        parameters["monoCriteriaInpuctDictionary"][monoCriterion]["criteria"],
                        parameters["monoCriteriaInpuctDictionary"][monoCriterion]["tolerance"],
                    )
                )
                criteriaInputs.append(
                    [parameters["monoCriteriaInpuctDictionary"][monoCriterion]["input"]]
                )
            # MultiCriterion
            multiCriterionInputDictionary = parameters["multiCriterionInputDictionary"]
            multiCriterionInputDictionary["criteria"] = criteriaArray
            multiCriterionInputDictionary["inputsForCriterion"] = criteriaInputs
            criterion = xmc.multiCriterion.MultiCriterion(**multiCriterionInputDictionary)
            # ErrorEstimator
            MSEErrorEstimator = xmc.errorEstimator.ErrorEstimator(
                **parameters["errorEstimatorInputDictionary"]
            )
            # HierarchyOptimiser
            hierarchyCostOptimiser = xmc.hierarchyOptimiser.HierarchyOptimiser(
                **parameters["hierarchyOptimiserInputDictionary"]
            )
            # EstimationAssembler
            if (
                "expectationAssembler"
                in parameters["estimationAssemblerInputDictionary"].keys()
            ):
                expectationAssembler = xmc.estimationAssembler.EstimationAssembler(
                    **parameters["estimationAssemblerInputDictionary"]["expectationAssembler"]
                )
            if (
                "discretizationErrorAssembler"
                in parameters["estimationAssemblerInputDictionary"].keys()
            ):
                discretizationErrorAssembler = xmc.estimationAssembler.EstimationAssembler(
                    **parameters["estimationAssemblerInputDictionary"][
                        "discretizationErrorAssembler"
                    ]
                )
            if "varianceAssembler" in parameters["estimationAssemblerInputDictionary"].keys():
                varianceAssembler = xmc.estimationAssembler.EstimationAssembler(
                    **parameters["estimationAssemblerInputDictionary"]["varianceAssembler"]
                )
            # MonteCarloSampler
            monteCarloSamplerInputDictionary = parameters["monteCarloSamplerInputDictionary"]
            monteCarloSamplerInputDictionary[
                "indexConstructorDictionary"
            ] = monteCarloIndexInputDictionary
            monteCarloSamplerInputDictionary["assemblers"] = [
                expectationAssembler,
                discretizationErrorAssembler,
                varianceAssembler,
            ]
            monteCarloSamplerInputDictionary["errorEstimators"] = [MSEErrorEstimator]
            mcSampler = xmc.monteCarloSampler.MonteCarloSampler(
                **monteCarloSamplerInputDictionary
            )
            # XMCAlgorithm
            XMCAlgorithmInputDictionary = parameters["XMCAlgorithmInputDictionary"]
            XMCAlgorithmInputDictionary["monteCarloSampler"] = mcSampler
            XMCAlgorithmInputDictionary["hierarchyOptimiser"] = hierarchyCostOptimiser
            XMCAlgorithmInputDictionary["stoppingCriterion"] = criterion
            algo = xmc.XMCAlgorithm(**XMCAlgorithmInputDictionary)

            if parameters["solverWrapperInputDictionary"]["asynchronous"] is True:
                algo.runAsynchronousXMC()
            else:
                algo.runXMC()

            # test
            estimations = get_value_from_remote(algo.estimation())
            estimated_mean = 1.47
            self.assertAlmostEqual(estimations[0], estimated_mean, delta=1.0)
            for level in algo.hierarchy():
                self.assertEqual(level[1], 15)

if __name__ == "__main__":
    unittest.main()
