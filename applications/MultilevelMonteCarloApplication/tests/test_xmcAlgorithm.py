"""
This test is designed to ensure the workflow Kratos-XMC works correctly under the following scenarios:
* workflow is serial,
* workflow is serial and managed by distributed environment scheduler.
To run the first scenario:
python3 test_xmcAlgorithm.py
To run with runcompss the second scenario:
sh test_runcompss_xmcALgorithm.sh
In this last case, the environment variable EXAQUTE_BACKEND has to be changed to pycompss; see the documentation related to the configuration of COMPSs for details.

Dependencies
------------
- KratosMultiphysics ≥ 9.0."Dev"-96fb824069, and applications:
   - ConvectionDiffusionApplication,
   - LinearSolversApplication,
   - MeshingApplication and
   - MultilevelMonteCarloApplication.
- COMPSs ≥ 2.8 (to run in parallel).
"""

# Import Python libraries
import unittest
import json
import sys
import os

# Import XMC, distributed environment
import xmc
from exaqute import get_value_from_remote


def isKratosFound():
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


def isMmgFound():
    try:
        import KratosMultiphysics
        import KratosMultiphysics.MeshingApplication

        if hasattr(KratosMultiphysics.MeshingApplication, "MmgProcess2D"):
            return True
    except ImportError:
        return False


class TestXMCAlgorithm(unittest.TestCase):
    def test_mc_Kratos(self):
        if not isKratosFound():
            self.skipTest(
                "Missing dependency: KratosMultiphysics or one of required applications. Check the test docstrings for details."
            )

        # read parameters
        parametersList = [
            "poisson_square_2d/problem_settings/parameters_xmc_test_mc_Kratos_asynchronous_poisson_2d.json",
            "poisson_square_2d/problem_settings/parameters_xmc_test_mc_Kratos_asynchronous_poisson_2d_RFF.json",
            "poisson_square_2d/problem_settings/parameters_xmc_test_mc_Kratos_asynchronous_poisson_2d_with_combined_power_sums.json",
            "poisson_square_2d/problem_settings/parameters_xmc_test_mc_Kratos_asynchronous_poisson_2d_with_10_combined_power_sums.json",
            "poisson_square_2d/problem_settings/parameters_xmc_test_mc_Kratos_poisson_2d.json",
            "poisson_square_2d/problem_settings/parameters_xmc_test_mc_Kratos_poisson_2d_with_combined_power_sums.json",
            "poisson_square_2d/problem_settings/poisson_multi-moment_mc.json",
            "poisson_square_2d/problem_settings/parameters_xmc_test_mc_Kratos_asynchronous_poisson_2d_fixedsamples.json",
            "poisson_square_2d/problem_settings/parameters_xmc_test_mc_Kratos_poisson_2d_fixedsamples.json",
        ]

        for parametersPath in parametersList:
            with open(parametersPath, "r") as parameter_file:
                parameters = json.load(parameter_file)
            # SolverWrapper
            parameters["solverWrapperInputDictionary"]["qoiEstimator"] = parameters[
                "monteCarloIndexInputDictionary"
            ]["qoiEstimator"]
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
            # MonoCriterion
            criteriaArray = []
            criteriaInputs = []
            for monoCriterion in parameters["monoCriteriaInputDictionary"]:
                criteriaArray.append(
                    xmc.monoCriterion.MonoCriterion(
                        parameters["monoCriteriaInputDictionary"][monoCriterion]["criteria"],
                        parameters["monoCriteriaInputDictionary"][monoCriterion]["tolerance"],
                    )
                )
                criteriaInputs.append(
                    [parameters["monoCriteriaInputDictionary"][monoCriterion]["input"]]
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
            if parameters["solverWrapperInputDictionary"]["asynchronous"]:
                self.assertEqual(algo.monteCarloSampler.samplesCounter,algo.hierarchy()[0][1])
                if parameters["samplerInputDictionary"]["randomGenerator"] == "xmc.randomGeneratorWrapper.EventDatabase":
                    self.assertEqual(algo.monteCarloSampler.batchIndices[-1][-1].sampler.randomGenerator._eventCounter,algo.hierarchy()[0][1])

    def test_mlmc_Kratos(self):
        if not isKratosFound():
            self.skipTest("Missing dependency: KratosMultiphysics or one of its applications")
        if not isMmgFound():
            self.skipTest(
                "Missing dependency: KratosMultiphysics.MeshingApplication with MMG support"
            )

        # read parameters
        parametersList = [
            "poisson_square_2d/problem_settings/parameters_xmc_test_mlmc_Kratos_asynchronous_poisson_2d.json",
            "poisson_square_2d/problem_settings/parameters_xmc_test_mlmc_Kratos_asynchronous_poisson_2d_with_combined_power_sums.json",
            "poisson_square_2d/problem_settings/parameters_xmc_test_mlmc_Kratos_poisson_2d.json",
            "poisson_square_2d/problem_settings/parameters_xmc_test_mlmc_Kratos_poisson_2d_with_combined_power_sums.json",
            "poisson_square_2d/problem_settings/poisson_multi-moment_mlmc.json",
            "poisson_square_2d/problem_settings/parameters_xmc_test_mlmc_Kratos_asynchronous_poisson_2d_with_combined_power_sums_multi.json",
            "poisson_square_2d/problem_settings/parameters_xmc_test_mlmc_Kratos_asynchronous_poisson_2d_with_combined_power_sums_multi_ensemble.json",
            "poisson_square_2d/problem_settings/parameters_xmc_test_mlmc_Kratos_asynchronous_poisson_2d_DAR.json",
            "poisson_square_2d/problem_settings/parameters_xmc_test_mlmc_Kratos_asynchronous_poisson_2d_fixedsamples.json",
        ]
        for parametersPath in parametersList:
            with open(parametersPath, "r") as parameter_file:
                parameters = json.load(parameter_file)
            # SolverWrapper
            parameters["solverWrapperInputDictionary"]["qoiEstimator"] = parameters[
                "monteCarloIndexInputDictionary"
            ]["qoiEstimator"]
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
            # MonoCriterion
            criteriaArray = []
            criteriaInputs = []
            for monoCriterion in parameters["monoCriteriaInputDictionary"]:
                criteriaArray.append(
                    xmc.monoCriterion.MonoCriterion(
                        parameters["monoCriteriaInputDictionary"][monoCriterion]["criteria"],
                        parameters["monoCriteriaInputDictionary"][monoCriterion]["tolerance"],
                    )
                )
                criteriaInputs.append(
                    [parameters["monoCriteriaInputDictionary"][monoCriterion]["input"]]
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
    unittest.main(verbosity=2)
