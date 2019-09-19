import KratosMultiphysics
import KratosMultiphysics.DEMApplication as DEMApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_guis
import test_particle_creator_destructor
import test_wall_creator_destructor
import test_analytics
import test_glued_particles
import test_restart
import test_DEM_2D

def AssembleTestSuites():

    ''' Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']
    smallSuite.addTest(test_guis.TestGUIs("test_GUIs_1"))
    smallSuite.addTest(test_guis.TestGUIs("test_GUIs_2"))
    smallSuite.addTest(test_particle_creator_destructor.TestParticleCreatorDestructor("test_CreateSphericParticle1"))
    smallSuite.addTest(test_particle_creator_destructor.TestParticleCreatorDestructor("test_CreateSphericParticle2"))
    smallSuite.addTest(test_wall_creator_destructor.TestWallCreatorDestructor("test_CreateWallTriangle"))
    smallSuite.addTest(test_analytics.TestAnalytics("test_Analytics_1"))
    #smallSuite.addTest(test_analytics.TestAnalytics("test_Analytics_2"))
    #smallSuite.addTest(test_analytics.TestAnalytics("test_Analytics_3"))
    smallSuite.addTest(test_glued_particles.TestGluedParticles("test_Glued_Particles_1"))
    smallSuite.addTest(test_DEM_2D.TestDEM2D("test_DEM2D_1"))
    smallSuite.addTest(test_restart.TestRestartOneBall("test_execution"))
    smallSuite.addTest(test_restart.TestRestartTwoBalls("test_execution"))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    # For very long tests that should not be in nightly and you can use to validate
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
