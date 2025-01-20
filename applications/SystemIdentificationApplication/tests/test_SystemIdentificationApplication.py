import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
import test_adjoint_sensors
import test_sensor_output_process
import test_system_identification
import responses.test_damage_response
import responses.test_sensor_count_response
import test_smooth_clamper
import test_mask_utils
import test_sensor_utils
import test_distance_matrix

def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites

    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_adjoint_sensors.TestDisplacementSensor]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_adjoint_sensors.TestStrainSensorShell]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_adjoint_sensors.TestStrainSensorSolids]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sensor_output_process.TestSensorOutputProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_smooth_clamper.TestSmoothClamper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mask_utils.TestMaskUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sensor_utils.TestSensorUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_distance_matrix.TestDistanceMatrix]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_system_identification.TestSystemIdentification]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_damage_response.TestDamageDetectionAdjointResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_damage_response.TestDamageDetectionResponse]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_sensor_count_response.TestSensorCountResponse]))

    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
