import KratosMultiphysics 
import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_input_output
import test_vms_adjoint_element_2d
import test_vms_sensitivity_for_one_time_step_2d
import test_vms_sensitivity_2d

## NIGTHLY TESTS

## VALIDATION TESTS

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

    # Create a test suite with the selected tests (Small tests):
    smallSuite = suites['small']
    smallSuite.addTest(test_input_output.TestCase('test_Execution'))
    smallSuite.addTest(test_vms_adjoint_element_2d.TestCase('test_MASS_MATRIX_0'))
    smallSuite.addTest(test_vms_adjoint_element_2d.TestCase('test_MASS_MATRIX_1'))
    smallSuite.addTest(test_vms_adjoint_element_2d.TestCase('test_ADJOINT_MATRIX_1'))
    smallSuite.addTest(test_vms_adjoint_element_2d.TestCase('test_ADJOINT_MATRIX_2'))
    smallSuite.addTest(test_vms_adjoint_element_2d.TestCase('test_SHAPE_DERIVATIVE_MATRIX_1'))
    smallSuite.addTest(test_vms_adjoint_element_2d.TestCase('test_SHAPE_DERIVATIVE_MATRIX_2'))
    smallSuite.addTest(test_vms_sensitivity_for_one_time_step_2d.TestCase('test_PrimalGradient'))
    smallSuite.addTest(test_vms_sensitivity_for_one_time_step_2d.TestCase('test_ElementSensitivity'))
    smallSuite.addTest(test_vms_sensitivity_for_one_time_step_2d.TestCase('test_AdjointBossakDragScheme'))
    smallSuite.addTest(test_vms_sensitivity_for_one_time_step_2d.TestCase('test_Sensitivity'))
    smallSuite.addTest(test_vms_sensitivity_2d.TestCase('test_OneElement'))
    smallSuite.addTest(test_vms_sensitivity_2d.TestCase('test_Cylinder'))


    # Create a test suite with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    
    # For very long tests that should not be in nightly and you can use to validate 
    validationSuite = suites['validation']

    # Create a test suite that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests([smallSuite])
    
    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
