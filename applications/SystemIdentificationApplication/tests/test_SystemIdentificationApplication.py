import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
import test_adjoint_sensors
import test_sensor_output_process
import test_system_identification
import responses.test_damage_response
import responses.test_temperature_response
import controls.test_data_values_control
import test_smooth_clamper
import test_mask_utils
import test_distance_matrix
import test_sensor_utils

def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites

    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_adjoint_sensors.TestDisplacementSensor]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_adjoint_sensors.TestStrainSensorShell]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_adjoint_sensors.TestStrainSensorSolids]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sensor_output_process.TestSensorOutputProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_smooth_clamper.TestSmoothClamper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_mask_utils.TestMaskUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_distance_matrix.TestDistanceMatrix]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sensor_utils.TestSensorUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_system_identification.TestSystemIdentification]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_damage_response.TestDamageDetectionAdjointResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_damage_response.TestDamageDetectionResponse]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_damage_response.TestDamageDetectionResponseStrainSensor]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_temperature_response.TestTemperatureDetectionAdjointResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_temperature_response.TestTemperatureDetectionResponse]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([responses.test_temperature_response.TestTemperatureDetectionResponseStrainSensor]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([controls.test_data_values_control.TestDataValuesControl]))

    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
