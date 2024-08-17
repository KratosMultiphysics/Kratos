import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
import test_adjoint_sensors
import test_sensor_output_process
import test_system_identification
import test_smooth_clamper
import test_distance_matrix
import test_mask_utils
import test_sensor_mask_status
import test_sensor_mask_status_kd_tree
import responses.test_sensor_count_response
import responses.test_sensor_distance_summation_response
import responses.test_sensor_isolation_response
import responses.test_sensor_coverage_response
import responses.test_sensor_localization_response

def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites

    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_adjoint_sensors.TestDisplacementSensor]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_adjoint_sensors.TestStrainSensor]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sensor_output_process.TestSensorOutputProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_smooth_clamper.TestSmoothClamper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_distance_matrix.TestDistanceMatrix]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mask_utils.TestMaskUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sensor_mask_status.TestSensorMaskStatus]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sensor_mask_status_kd_tree.TestSensorMaskStatusKDTree]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_system_identification.TestSystemIdentification]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_sensor_count_response.TestSensorCountResponse]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_sensor_distance_summation_response.TestSensorDistanceSummationResponse]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_sensor_isolation_response.TestSensorIsolationResponse]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_sensor_coverage_response.TestSensorCoverageResponse]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_sensor_localization_response.TestSensorLocalizationResponse]))

    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
