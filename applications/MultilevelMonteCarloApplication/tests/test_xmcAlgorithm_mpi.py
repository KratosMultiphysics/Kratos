"""
This test is designed to ensure the workflow Kratos-XMC works correctly under the following scenarios:
* workflow is MPI parallel,
* workflow is serial, only Kratos tasks are MPI parallel and are scheduled by distributed environment scheduler.
To run the first scenario:
mpirun -n $number_processes python3 test_xmcAlgorithm_mpi.py
To run with runcompss the second scenario:
sh test_runcompss_xmcALgorithm_mpi.sh
In the latter case, the environment variable EXAQUTE_BACKEND has to be set to pycompss,
and in test_runcompss_xmcALgorithm_mpi.sh the test to run must be selected.
In the former case, the environment variable EXAQUTE_BACKEND has to be set to local.
See the documentation section related to the configuration of COMPSs for details.

Dependencies
------------
- KratosMultiphysics ≥ 9.0."Dev"-96fb824069 configured with the CMAKE flag USE_MPI set to "ON",
  and applications:
   - FluidDynamicsApplication,
   - LinearSolversApplication,
   - MappingApplication,
   - MeshingApplication,
   - MetisApplication,
   - MultilevelMonteCarloApplication,
   - StatisticsApplication,
   - TrilinosApplication.
- COMPSs ≥ 2.8 (to run in parallel) with its optional dependency: mpi4py.
- ParMmg ≥ commit 5ffc6ada4afb1af50a43e1fa6f4c409cff2ea25c (as of 2020-02-16, this is branch 'develop' of https://github.com/MmgTools/parmmg) for "mpi_test_mlmc_Kratos_ParMmg" test.
  See https://github.com/KratosMultiphysics/Kratos/wiki/%5BUtilities%5D-ParMmg-Process
  for compiling KratosMultiphysics with ParMmg support.
"""

# Import Python libraries
import unittest
import json
import sys
import os

# Import XMC, distributed environment API
import xmc
from exaqute import get_value_from_remote
import xmc.methodDefs_momentEstimator.computeCentralMoments as ccm


class WorkFolderScope:
    """
    This class is taken from the core of Kratos Multiphysics.

    Helper-class to execute test in a specific target path
    Input:
    - rel_path_work_folder: String
        Relative path of the target dir from the calling script

    - file_path: String
        Absolute path of the calling script

    - add_to_path: Bool
        "False" (default) if no need to add the target dir to the path, "True" otherwise.
    """

    def __init__(self, rel_path_work_folder, file_path, add_to_path=False):
        self.currentPath = os.getcwd()
        self.add_to_path = add_to_path
        if self.add_to_path:
            self.currentPythonpath = sys.path
        self.scope = os.path.abspath(
            os.path.join(os.path.dirname(os.path.realpath(file_path)), rel_path_work_folder)
        )

    def __enter__(self):
        os.chdir(self.scope)
        if self.add_to_path:
            sys.path.append(self.scope)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.currentPath)
        if self.add_to_path:
            sys.path.remove(self.scope)


try:
    from KratosMultiphysics.kratos_utilities import (
        CheckIfApplicationsAvailable,
        IsMPIAvailable,
    )

    is_Kratos = (
        CheckIfApplicationsAvailable(
            "FluidDynamicsApplication",
            "LinearSolversApplication",
            "MappingApplication",
            "MeshingApplication",
            "MetisApplication",
            "MultilevelMonteCarloApplication",
            "StatisticsApplication",
            "TrilinosApplication",
        )
        and IsMPIAvailable()
    )
except ImportError:
    is_Kratos = False

try:
    import KratosMultiphysics
    import KratosMultiphysics.MeshingApplication

    is_ParMmg = hasattr(KratosMultiphysics.MeshingApplication, "ParMmgProcess3D")
except ImportError:
    is_ParMmg = False


class TestXMCAlgorithmMPI(unittest.TestCase):
    @unittest.skipIf(
        not is_Kratos,
        "Missing dependency: KratosMultiphysics, MPI or one of required applications. Check the test docstrings for details.",
    )
    def mpi_test_mc_Kratos(self):

        # read parameters
        parametersList = [
            "problem_settings/parameters_xmc_asynchronous_mc_SARTAAO.json",
            "problem_settings/parameters_xmc_asynchronous_mc_SARMT.json",
            "problem_settings/parameters_xmc_asynchronous_mc_DAR.json",
            "problem_settings/parameters_xmc_asynchronous_mc_RFF.json",
        ]

        with WorkFolderScope("caarc_wind_mpi/", __file__, add_to_path=True):
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
                monteCarloIndexInputDictionary[
                    "samplerInputDictionary"
                ] = samplerInputDictionary
                # MonoCriterion
                criteriaArray = []
                criteriaInputs = []
                for monoCriterion in parameters["monoCriteriaInputDictionary"]:
                    criteriaArray.append(
                        xmc.monoCriterion.MonoCriterion(
                            parameters["monoCriteriaInputDictionary"][monoCriterion][
                                "criteria"
                            ],
                            parameters["monoCriteriaInputDictionary"][monoCriterion][
                                "tolerance"
                            ],
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
                        **parameters["estimationAssemblerInputDictionary"][
                            "expectationAssembler"
                        ]
                    )
                if (
                    "varianceAssembler"
                    in parameters["estimationAssemblerInputDictionary"].keys()
                ):
                    varianceAssembler = xmc.estimationAssembler.EstimationAssembler(
                        **parameters["estimationAssemblerInputDictionary"]["varianceAssembler"]
                    )
                # MonteCarloSampler
                monteCarloSamplerInputDictionary = parameters[
                    "monteCarloSamplerInputDictionary"
                ]
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
                self.assertGreater(-estimations[0], 0)
                self.assertEqual(algo.hierarchy()[0][1], 10)
                # check moment estimator
                sample_counter = (
                    algo.monteCarloSampler.indices[0].qoiEstimator[0]._sampleCounter
                )
                S1 = get_value_from_remote(
                    algo.monteCarloSampler.indices[0].qoiEstimator[0].powerSums[0][0]
                )
                h1 = get_value_from_remote(
                    ccm.computeCentralMomentsOrderOneDimensionZero(S1, sample_counter)
                )
                self.assertGreater(-h1, 0)
                self.assertEqual(sample_counter, 10)
                # check multi moment estimator
                sample_counter = (
                    algo.monteCarloSampler.indices[0].qoiEstimator[1]._sampleCounter
                )
                self.assertEqual(sample_counter, 10)
                # check combined moment estimator
                sample_counter = (
                    algo.monteCarloSampler.indices[0].qoiEstimator[2]._sampleCounter
                )
                S1 = get_value_from_remote(
                    algo.monteCarloSampler.indices[0].qoiEstimator[2].powerSums[0][0]
                )
                h1 = get_value_from_remote(
                    ccm.computeCentralMomentsOrderOneDimensionZero(S1, sample_counter)
                )
                self.assertGreater(-h1, 0)
                self.assertEqual(sample_counter, 10)
                # check multi combined moment estimator
                sample_counter = (
                    algo.monteCarloSampler.indices[0].qoiEstimator[3]._sampleCounter
                )
                self.assertEqual(sample_counter, 10)

    @unittest.skipIf(
        not is_Kratos,
        "Missing dependency: KratosMultiphysics, MPI or one of required applications. Check the test docstrings for details.",
    )
    def mpi_test_mlmc_Kratos(self):

        # read parameters
        parametersList = ["problem_settings/parameters_xmc_asynchronous_mlmc_RFF.json"]

        with WorkFolderScope("caarc_wind_mpi/", __file__, add_to_path=True):
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
                monteCarloIndexInputDictionary[
                    "samplerInputDictionary"
                ] = samplerInputDictionary
                # MonoCriterion
                criteriaArray = []
                criteriaInputs = []
                for monoCriterion in parameters["monoCriteriaInputDictionary"]:
                    criteriaArray.append(
                        xmc.monoCriterion.MonoCriterion(
                            parameters["monoCriteriaInputDictionary"][monoCriterion][
                                "criteria"
                            ],
                            parameters["monoCriteriaInputDictionary"][monoCriterion][
                                "tolerance"
                            ],
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
                monteCarloSamplerInputDictionary = parameters[
                    "monteCarloSamplerInputDictionary"
                ]
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
                estimations = get_value_from_remote(algo.estimation())
                self.assertGreater(sum(estimations), 0)
                for level in algo.hierarchy():
                    self.assertEqual(level[1], 10)
                # check moment estimator - level 0
                sample_counter = (
                    algo.monteCarloSampler.indices[0].qoiEstimator[0]._sampleCounter
                )
                S1 = get_value_from_remote(
                    algo.monteCarloSampler.indices[0].qoiEstimator[0].powerSums[0][0]
                )
                h1 = get_value_from_remote(
                    ccm.computeCentralMomentsOrderOneDimensionZero(S1, sample_counter)
                )
                self.assertGreater(-h1, 0)
                self.assertEqual(sample_counter, 10)
                # check multi moment estimator - level 1
                sample_counter = (
                    algo.monteCarloSampler.indices[1].qoiEstimator[1]._sampleCounter
                )
                self.assertEqual(sample_counter, 10)
                # check combined moment estimator - level 2
                sample_counter = (
                    algo.monteCarloSampler.indices[2].qoiEstimator[2]._sampleCounter
                )
                S1 = get_value_from_remote(
                    algo.monteCarloSampler.indices[2].qoiEstimator[2].powerSums[0][0]
                )
                h1 = get_value_from_remote(
                    ccm.computeCentralMomentsOrderOneDimensionZero(S1, sample_counter)
                )
                self.assertEqual(sample_counter, 10)
                # check multi combined moment estimator - level 2
                sample_counter = (
                    algo.monteCarloSampler.indices[2].qoiEstimator[3]._sampleCounter
                )
                self.assertEqual(sample_counter, 10)

    @unittest.skipIf(
        not is_Kratos,
        "Missing dependency: KratosMultiphysics, MPI or one of required applications. Check the test docstrings for details.",
    )
    @unittest.skipIf(
        not is_ParMmg,
        "Missing dependency: ParMmg. You need to compile Kratos with ParMmg to run this test.",
    )
    def mpi_test_mlmc_Kratos_ParMmg(self):

        # read parameters
        parametersList = [
            "problem_settings/parameters_xmc_asynchronous_mlmc_SAR.json",
            "problem_settings/parameters_xmc_asynchronous_mlmc_DAR.json",
        ]

        with WorkFolderScope("caarc_wind_mpi/", __file__, add_to_path=True):
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
                monteCarloIndexInputDictionary[
                    "samplerInputDictionary"
                ] = samplerInputDictionary
                # MonoCriterion
                criteriaArray = []
                criteriaInputs = []
                for monoCriterion in parameters["monoCriteriaInputDictionary"]:
                    criteriaArray.append(
                        xmc.monoCriterion.MonoCriterion(
                            parameters["monoCriteriaInputDictionary"][monoCriterion][
                                "criteria"
                            ],
                            parameters["monoCriteriaInputDictionary"][monoCriterion][
                                "tolerance"
                            ],
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
                monteCarloSamplerInputDictionary = parameters[
                    "monteCarloSamplerInputDictionary"
                ]
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
                estimations = get_value_from_remote(algo.estimation())
                self.assertGreater(sum(estimations), 0)
                for level in algo.hierarchy():
                    self.assertEqual(level[1], 5)
                # check moment estimator - level 0
                sample_counter = (
                    algo.monteCarloSampler.indices[0].qoiEstimator[0]._sampleCounter
                )
                S1 = get_value_from_remote(
                    algo.monteCarloSampler.indices[0].qoiEstimator[0].powerSums[0][0]
                )
                h1 = get_value_from_remote(
                    ccm.computeCentralMomentsOrderOneDimensionZero(S1, sample_counter)
                )
                self.assertGreater(-h1, 0)
                self.assertEqual(sample_counter, 5)
                # check multi moment estimator - level 1
                sample_counter = (
                    algo.monteCarloSampler.indices[1].qoiEstimator[1]._sampleCounter
                )
                self.assertEqual(sample_counter, 5)
                # check combined moment estimator - level 2
                sample_counter = (
                    algo.monteCarloSampler.indices[2].qoiEstimator[2]._sampleCounter
                )
                S1 = get_value_from_remote(
                    algo.monteCarloSampler.indices[2].qoiEstimator[2].powerSums[0][0]
                )
                h1 = get_value_from_remote(
                    ccm.computeCentralMomentsOrderOneDimensionZero(S1, sample_counter)
                )
                self.assertEqual(sample_counter, 5)
                # check multi combined moment estimator - level 2
                sample_counter = (
                    algo.monteCarloSampler.indices[2].qoiEstimator[3]._sampleCounter
                )
                self.assertEqual(sample_counter, 5)


if __name__ == "__main__":
    # Define a loader to catch all and only MPI tests
    mpi_loader = unittest.TestLoader()
    mpi_loader.testMethodPrefix = "mpi_test_"
    unittest.main(testLoader=mpi_loader, verbosity=2)
