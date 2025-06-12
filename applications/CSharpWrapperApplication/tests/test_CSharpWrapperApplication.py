# import Kratos
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.CSharpWrapperApplication as CSharpWrapperApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import subprocess
import subprocess

# Using kratos_utilities
import KratosMultiphysics.kratos_utilities as kratos_utilities
if kratos_utilities.CheckIfApplicationsAvailable("ExternalSolversApplication"):
    has_external_solvers_application = True
else:
    has_external_solvers_application = False

# Import the tests o test_classes to create the suits
## SMALL TESTS

## NIGTHLY TESTS

## VALIDATION TESTS

def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']

    # Create a test suit with the selected tests plus all small tests
    nightlySuite = suites['nightly']

    ### BEGIN SMALL SUITE ###

    ### END SMALL SUITE ###

    ### BEGIN NIGHTLY SUITE ###

    ### END VALIDATION SUITE ###

    ### BEGIN VALIDATION SUITE ###

    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTests(nightlySuite)

    ### END VALIDATION ###

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightlySuite) # Already contains the smallSuite
    validationSuite.addTests(allSuite) # Validation contains all

    # Manual list for debugging
    #allSuite.addTests(
        #KratosUnittest.TestLoader().loadTestsFromTestCases([
            #### STANDALONE
            #### SMALL
            #### NIGTHLY
            #### VALIDATION
        #])
    #)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
