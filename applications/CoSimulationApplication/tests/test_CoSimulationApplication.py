import KratosMultiphysics.KratosUnittest as KratosUnittest

from convergence_accelerators.test_iqni import TestConvergenceAcceleratorIQNI
from convergence_criteria.test_absolute_norm import TestConvergenceCriterionAbsoluteNorm
from convergence_criteria.test_and import TestConvergenceCriterionAnd
from convergence_criteria.test_iteration_limit import TestConvergenceCriterionIterationLimit
from convergence_criteria.test_or import TestConvergenceCriterionOr
from convergence_criteria.test_relative_norm import TestConvergenceCriterionRelativeNorm
from mappers.test_linear import TestMapperLinear
from predictors.test_linear import TestPredictorLinear
from solver_wrappers.pipe.test_flow_solver import TestSolverWrapperPipeFlowSolver
from solver_wrappers.pipe.test_structure_solver import TestSolverWrapperPipeStructureSolver
from pykratos.test_parameters_class import TestPyKratosParameters
from pykratos.test_variables import TestPyKratosVariables
from pykratos.test_cosimulation_interface import TestCoSimulationInterface


def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    suites = KratosUnittest.KratosSuites

    smallSuite = suites['small']  # These tests are executed by the continuous integration tool
    smallSuite.addTest(TestConvergenceAcceleratorIQNI("test_convergence_accelerator_iqni"))
    smallSuite.addTest(TestConvergenceCriterionAbsoluteNorm("test_convergence_criterion_absolute_norm"))
    smallSuite.addTest(TestConvergenceCriterionAnd("test_convergence_criterion_and"))
    smallSuite.addTest(TestConvergenceCriterionIterationLimit("test_convergence_criterion_iteration_limit"))
    smallSuite.addTest(TestConvergenceCriterionOr("test_convergence_criterion_or"))
    smallSuite.addTest(TestConvergenceCriterionRelativeNorm("test_convergence_criterion_relative_norm"))
    smallSuite.addTest(TestMapperLinear("test_mapper_linear"))
    smallSuite.addTest(TestPredictorLinear("test_predictor_linear"))
    smallSuite.addTest(TestSolverWrapperPipeFlowSolver("test_solver_wrapper_pipe_flow_solver"))
    smallSuite.addTest(TestSolverWrapperPipeStructureSolver("test_solver_wrapper_pipe_structure_solver"))
    smallSuite.addTest(TestPyKratosParameters("test_pykratos_parameters"))
    smallSuite.addTest(TestPyKratosVariables("test_pykratos_variables"))
    smallSuite.addTest(TestCoSimulationInterface("test_cosimulation_interface"))

    nightlySuite = suites['nightly']  # These tests are executed in the nightly build
    nightlySuite.addTests(smallSuite)

    validationSuite = suites['validation']   # These tests are very long and should not be in nightly, for validation

    # Create a test suit that contains all tests:
    allSuite = suites['all']
    allSuite.addTests(nightlySuite)  # Already contains the smallSuite
    allSuite.addTests(validationSuite)

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
