# Definition of the classes for the NIGHTLY TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class candelier_no_history_test(TF.TestFactory):
     file_name = "candelier_tests/candelier"
     file_parameters = "candelier_tests/ProjectParametersDEM.json"
     file_fluid_parameters = "candelier_tests/ProjectParametersFluid.json"
     file_dem_parameters = "candelier_tests/ProjectParametersDEMDEM.json"

def SetTestSuite(suites):
    night_suite = suites['nightly']

    night_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            candelier_no_history_test
            ])
    )

    return night_suite
