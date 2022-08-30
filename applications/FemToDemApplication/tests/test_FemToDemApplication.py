# import Kratos
import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Small femdem cases
import main_coupling_for_testing
import main_coupling_total_lagrangian_for_testing
import main_coupling_for_testing_face_load
import main_coupling_for_testing_tables
import main_coupling_fracture_3_point_for_testing
import main_coupling_free_fall_testing
import main_fsi_aitken_for_testing

def AssembleTestSuites():

    suites = KratosUnittest.KratosSuites

    smallSuite = suites['small']

    smallSuite.addTest(main_coupling_for_testing.TestAnalytics("small_strain")) #defined inside main_coupling_for_testing
    smallSuite.addTest(main_coupling_total_lagrangian_for_testing.TestAnalytics("total_lagrangian")) #defined inside main_coupling_total_lagrangian_for_testing
    smallSuite.addTest(main_coupling_for_testing_tables.TestAnalytics("tables")) #defined inside main_coupling_for_testing_tables
    smallSuite.addTest(main_coupling_fracture_3_point_for_testing.TestAnalytics("fracture_3_point")) #defined inside main_coupling_fracture_3_point_for_testing
    smallSuite.addTest(main_coupling_free_fall_testing.TestAnalytics("free_fall")) #defined inside main_coupling_fracture_3_point_for_testing
    smallSuite.addTest(main_fsi_aitken_for_testing.TestAnalytics("two_dimensional_fsi")) #defined inside main_fsi_aitken_for_testing
    smallSuite.addTest(main_coupling_for_testing_face_load.TestAnalytics("face_load")) #defined inside main_coupling_for_testing_face_load

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']

    # nightSuite.addTests(smallSuite)

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
