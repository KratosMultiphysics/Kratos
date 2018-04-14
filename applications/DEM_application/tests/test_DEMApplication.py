import KratosMultiphysics 
import KratosMultiphysics.DEMApplication as DEMApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_particle_creator_destructor


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
    smallSuite.addTest(test_particle_creator_destructor.TestParticleCreatorDestructor("test_CreateSphericParticle1"))
    smallSuite.addTest(test_particle_creator_destructor.TestParticleCreatorDestructor("test_CreateSphericParticle2"))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    
    # For very long tests that should not be in nighly and you can use to validate 
    validationSuite = suites['validation']
    
    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        smallSuite
        #KratosUnittest.TestLoader().loadTestsFromTestCases([])
    )
    

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
