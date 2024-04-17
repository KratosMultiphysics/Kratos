import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
import test_adjoint_sensors
import test_sensor_output_process
import test_system_identification
import test_smooth_clamper
import test_mask_utils
import responses.test_sensor_coverage_response
import responses.test_sensor_count_response
import responses.test_sensor_localization_response
import responses.test_sensor_inverse_distance_response
import test_sensor_data_model_part_controller
import controls.test_sensor_status_control

def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites

    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_adjoint_sensors.TestDisplacementSensor]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_adjoint_sensors.TestStrainSensor]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sensor_output_process.TestSensorOutputProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_smooth_clamper.TestSmoothClamper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mask_utils.TestMaskUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_sensor_coverage_response.TestSensorCoverageResponse]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_sensor_count_response.TestSensorCountResponse]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_sensor_localization_response.TestSensorLocalizationResponse]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_sensor_inverse_distance_response.TestSensorInverseDistanceResponse]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sensor_data_model_part_controller.TestSensorDataModelPartController]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([controls.test_sensor_status_control.TestSensorStatusControl]))
    # smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_system_identification.TestSystemIdentification]))

    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
