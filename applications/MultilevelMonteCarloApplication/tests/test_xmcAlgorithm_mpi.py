"""
This test is designed to ensure the workflow Kratos-XMC works correctly under the following scenarios:
* workflow is MPI parallel,
* workflow is serial, only Kratos tasks are MPI parallel and are scheduled by distributed environment scheduler.
To run the first scenario:
mpirun -n $number_processes python3 test_xmcAlgorithm_mpi.py
To run with runcompss the second scenario:
sh test_runcompss_xmcALgorithm_mpi.sh
In this last case, the appropriate import has to be changed in xmc/distributedEnvironmentFramework.py.

Required KratosMultiphysics application are: "ExaquteSandboxApplication", "FluidDynamicsApplication", "LinearSolversApplication", "MappingApplication", "MeshingApplication", "MetisApplication", "MultilevelMonteCarloApplication", "StatisticsApplication", "TrilinosApplication".
It is required to configure Kratos with the CMAKE flag USE_MPI set to "ON".
"""

# Import Python libraries
import unittest
import json
import sys
import os

# Import XMC, distributed environment API
import xmc
from xmc.distributedEnvironmentFramework import get_value_from_remote
import xmc.methodDefs_momentEstimator.computeCentralMoments as ccm

def isKratosFound():
    try:
        from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable, IsMPIAvailable
        return CheckIfApplicationsAvailable("ExaquteSandboxApplication", "FluidDynamicsApplication", "LinearSolversApplication", "MappingApplication", "MeshingApplication", "MetisApplication", "MultilevelMonteCarloApplication", "StatisticsApplication", "TrilinosApplication") and IsMPIAvailable()
    except ImportError:
        return False


class TestXMCAlgorithmMPI(unittest.TestCase):
    def test_mc_Kratos(self):
        if not isKratosFound():
            self.skipTest("Missing dependency: KratosMultiphysics, MPI or one of required applications. Check the test docstrings for details.")

        # read parameters
        parametersList = [
            "rectangle_wind_mpi/problem_settings/parameters_xmc_asynchronous_mc_SARTAAO.json",
            "rectangle_wind_mpi/problem_settings/parameters_xmc_asynchronous_mc_SARMT.json",
            "rectangle_wind_mpi/problem_settings/parameters_xmc_asynchronous_mc_DAR.json",
            "rectangle_wind_mpi/problem_settings/parameters_xmc_asynchronous_mc_RFF.json"
        ]

        for parametersPath in parametersList:
            with open(parametersPath, "r") as parameter_file:
                parameters = json.load(parameter_file)
            # add path of the problem folder to python path
            problem_id = parameters["solverWrapperInputDictionary"]["problemId"]
            sys.path.append(
                os.path.join(problem_id)
            )
            # SolverWrapper
            parameters["solverWrapperInputDictionary"]["qoiEstimator"] = parameters["monteCarloIndexInputDictionary"]["qoiEstimator"]
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
            # such results are not accurate, since we run the problem for few decimals
            # and coarse meshes instead of hundreds of seconds and finer meshes
            # check quantity of interest used to check convergence
            estimations = get_value_from_remote(algo.estimation())
            estimated_mean = 143
            self.assertAlmostEqual(estimations[0], estimated_mean, delta=10.0)
            self.assertEqual(algo.hierarchy()[0][1], 10)
            # check moment estimator
            sample_counter = algo.monteCarloSampler.indices[0].qoiEstimator[0]._sampleCounter
            S1 = get_value_from_remote(algo.monteCarloSampler.indices[0].qoiEstimator[0].powerSums[0][0])
            h1 = get_value_from_remote(ccm.computeCentralMomentsOrderOneDimensionZero(S1,sample_counter))
            self.assertAlmostEqual(h1, 143, delta=10.0)
            self.assertEqual(sample_counter, 10)
            # check multi moment estimator
            sample_counter = algo.monteCarloSampler.indices[0].qoiEstimator[1]._sampleCounter
            self.assertEqual(sample_counter, 10)
            # check combined moment estimator
            sample_counter = algo.monteCarloSampler.indices[0].qoiEstimator[2]._sampleCounter
            S1 = get_value_from_remote(algo.monteCarloSampler.indices[0].qoiEstimator[2].powerSums[0][0])
            h1 = get_value_from_remote(ccm.computeCentralMomentsOrderOneDimensionZero(S1,sample_counter))
            self.assertAlmostEqual(h1, 143, delta=10.0)
            self.assertEqual(sample_counter, 60)
            # check multi combined moment estimator
            sample_counter = algo.monteCarloSampler.indices[0].qoiEstimator[3]._sampleCounter
            self.assertEqual(sample_counter, 60)


    def test_mlmc_Kratos(self):
        if not isKratosFound():
            self.skipTest("Missing dependency: KratosMultiphysics, MPI or one of required applications. Check the test docstrings for details.")
        # here we should skip if ParMMg is not found, for SAR and DAR.

        # read parameters
        parametersList = [
            "rectangle_wind_mpi/problem_settings/parameters_xmc_asynchronous_mlmc_RFF.json"
        ]

        for parametersPath in parametersList:
            with open(parametersPath, "r") as parameter_file:
                parameters = json.load(parameter_file)
            # add path of the problem folder to python path
            problem_id = parameters["solverWrapperInputDictionary"]["problemId"]
            sys.path.append(
                os.path.join(problem_id)
            )
            # SolverWrapper
            parameters["solverWrapperInputDictionary"]["qoiEstimator"] = parameters["monteCarloIndexInputDictionary"]["qoiEstimator"]
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
            # MonteCarloSampler
            monteCarloSamplerInputDictionary = parameters["monteCarloSamplerInputDictionary"]
            monteCarloSamplerInputDictionary[
                "indexConstructorDictionary"
            ] = monteCarloIndexInputDictionary
            monteCarloSamplerInputDictionary["errorEstimators"] = [MSEErrorEstimator]
            # EstimationAssembler
            monteCarloSamplerInputDictionary["assemblers"] = []
            for key, dicArgs in parameters["estimationAssemblerInputDictionary"].items():
                monteCarloSamplerInputDictionary["assemblers"].append(
                    xmc.estimationAssembler.EstimationAssembler(**dicArgs)
                )
            mcSampler = xmc.monteCarloSampler.MonteCarloSampler(**monteCarloSamplerInputDictionary)
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
            # such results are not accurate, since we run the problem for few decimals
            # and coarse meshes instead of hundreds of seconds and finer meshes
            estimations = get_value_from_remote(algo.estimation())
            estimated_mean = 220
            self.assertAlmostEqual(sum(estimations), estimated_mean, delta=15.0)
            for level in algo.hierarchy():
                self.assertEqual(level[1], 10)
            # check moment estimator - level 0
            sample_counter = algo.monteCarloSampler.indices[0].qoiEstimator[0]._sampleCounter
            S1 = get_value_from_remote(algo.monteCarloSampler.indices[0].qoiEstimator[0].powerSums[0][0])
            h1 = get_value_from_remote(ccm.computeCentralMomentsOrderOneDimensionZero(S1,sample_counter))
            self.assertAlmostEqual(h1, 143, delta=15.0)
            self.assertEqual(sample_counter, 10)
            # check multi moment estimator - level 1
            sample_counter = algo.monteCarloSampler.indices[1].qoiEstimator[1]._sampleCounter
            self.assertEqual(sample_counter, 10)
            # check combined moment estimator - level 2
            sample_counter = algo.monteCarloSampler.indices[2].qoiEstimator[2]._sampleCounter
            S1 = get_value_from_remote(algo.monteCarloSampler.indices[2].qoiEstimator[2].powerSums[0][0])
            h1 = get_value_from_remote(ccm.computeCentralMomentsOrderOneDimensionZero(S1,sample_counter))
            self.assertAlmostEqual(h1, -950, delta=15.0)
            self.assertEqual(sample_counter, 60)
            # check multi combined moment estimator - level 2
            sample_counter = algo.monteCarloSampler.indices[2].qoiEstimator[3]._sampleCounter
            self.assertEqual(sample_counter, 60)

if __name__ == "__main__":
    unittest.main()
