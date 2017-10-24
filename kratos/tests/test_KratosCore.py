from __future__ import print_function, absolute_import, division

# import Kratos
from KratosMultiphysics import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suites
#from test_kratos_parameters import TestParameters as TParameters
from test_model_part_io import TestModelPartIO as TModelPartIO
from test_model_part import TestModelPart as TModelPart
import test_kratos_parameters
import test_materials_input
import test_geometries
import test_linear_solvers
import test_eigen_solvers
import test_condition_number
import test_processes
import test_importing
import test_connectivity_preserve_modeler
import test_model
import test_redistance
import test_variable_utils
import test_reorder
import test_exact_integration


def AssambleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    suites = KratosUnittest.KratosSuites

    # Create a test suite with the selected tests (Small tests):
    smallSuite = suites['small']
    #smallSuite.addTest(TModelPartIO('test_model_part_io_read_model_part'))
    #smallSuite.addTest(TModelPartIO('test_model_part_io_write_model_part'))
    #smallSuite.addTest(TModelPart('test_model_part_properties'))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TModelPart]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TModelPartIO]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_materials_input.TestMaterialsInput]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_geometries.TestGeometry]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_kratos_parameters.TestParameters]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_linear_solvers.TestLinearSolvers]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_eigen_solvers.TestEigenSolvers]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_condition_number.TestConditionNumber]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_processes.TestProcesses]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_importing.TestImporting]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_connectivity_preserve_modeler.TestConnectivityPreserveModeler]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_model.TestModel]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_redistance.TestRedistance]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_variable_utils.TestVariableUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_reorder.TestReorder]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_exact_integration.TestExactIntegration]))

    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']

    nightSuite.addTests(map(TModelPart, [
        'test_model_part_sub_model_parts',
        'test_model_part_nodes',
        'test_model_part_tables'
    ]))

    #nightSuite.addTests(map(TParameters, [
        #'test_kratos_parameters',
        #'test_kratos_change_parameters',
        #'test_kratos_copy_parameters',
        #'test_kratos_wrong_parameters'
    #]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TModelPartIO]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_kratos_parameters.TestParameters]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_materials_input.TestMaterialsInput]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_geometries.TestGeometry]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_linear_solvers.TestLinearSolvers]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_eigen_solvers.TestEigenSolvers]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_condition_number.TestConditionNumber]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_processes.TestProcesses]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_importing.TestImporting]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_connectivity_preserve_modeler.TestConnectivityPreserveModeler]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_model.TestModel]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_redistance.TestRedistance]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_variable_utils.TestVariableUtils]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_reorder.TestReorder]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_exact_integration.TestExactIntegration]))


    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            TModelPartIO,
            TModelPart,
            test_kratos_parameters.TestParameters,
            test_materials_input.TestMaterialsInput,
            test_geometries.TestGeometry,
            test_linear_solvers.TestLinearSolvers,
            test_eigen_solvers.TestEigenSolvers,
            test_condition_number.TestConditionNumber,
            test_processes.TestProcesses,
            test_importing.TestImporting,
            test_connectivity_preserve_modeler.TestConnectivityPreserveModeler,
            test_model.TestModel,
            test_redistance.TestRedistance,
            test_variable_utils.TestVariableUtils,
            test_reorder.TestReorder,
            test_exact_integration.TestExactIntegration
        ])
    )

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
