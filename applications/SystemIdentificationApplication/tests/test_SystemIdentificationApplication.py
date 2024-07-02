import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
import test_adjoint_sensors
import test_sensor_output_process
import test_system_identification
import test_smooth_clamper

def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites

    smallSuite = suites['small']
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_adjoint_sensors.TestDisplacementSensor]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_adjoint_sensors.TestStrainSensor]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_sensor_output_process.TestSensorOutputProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_smooth_clamper.TestSmoothClamper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([test_system_identification.TestSystemIdentification]))

    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)

    allSuite = suites['all']
    allSuite.addTests(nightSuite)

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
