# Import Python libraries
import unittest
import json
import sys
import os

# Import XMC, distributed environment API
import xmc
from xmc.distributedEnvironmentFramework import get_value_from_remote

def kratosFound():
    try:
        import KratosMultiphysics
        import KratosMultiphysics.mpi
        return True
    except ImportError:
        return False


class TestXMCAlgorithmMPI(unittest.TestCase):
    def test_mc_Kratos(self):
        if not kratosFound():
            self.skipTest("Missing dependencies: KratosMultiphysics and/or KratosMultiphysics.mpi")

        # read parameters
        parametersList = [
            "test_xmcAlgorithm_mpi/problem_settings/parameters_xmc_asynchronous_mc_SARTAAO.json",
            "test_xmcAlgorithm_mpi/problem_settings/parameters_xmc_asynchronous_mc_SARMT.json",
            "test_xmcAlgorithm_mpi/problem_settings/parameters_xmc_asynchronous_mc_DAR.json",
            "test_xmcAlgorithm_mpi/problem_settings/parameters_xmc_asynchronous_mc_RFF.json"
        ]

        for parametersPath in parametersList:
            with open(parametersPath, "r") as parameter_file:
                parameters = json.load(parameter_file)
            # add path of the problem folder to python path
            problem_id = parameters["solverWrapperInputDictionary"]["problemId"]
            sys.path.append(
                os.path.join(problem_id)
            )
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
            estimated_mean = 143
            self.assertAlmostEqual(estimations[0], estimated_mean, delta=10.0) # such result is not accurate, since we run the problem for few decimals instead of hundreds of seconds
            self.assertEqual(algo.hierarchy()[0][1], 10)

if __name__ == "__main__":
    unittest.main()
