# Definition of the classes for the VALIDATION TESTS

#Iimport Kratos
import KratosMultiphysics
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.SwimmingDEMApplication

# Import TestFactory
import FluidDEMTestFactory as FDEMTF
import SPFEMTestFactory as SPFEMTF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class fluid_dem_coupling_one_way_test(FDEMTF.TestFactory):
     file_name = "fluid_dem_tests/settling_cube"
     file_parameters = "fluid_dem_tests/ProjectParameters.json"

class sdem_pfem_coupling_one_way_test(SPFEMTF.TestFactory):
     file_name = "PFEM-DEM_tests/sdem_pfem_coupling_one_way_test"
     file_parameters = "PFEM-DEM_tests/ProjectParameters.json"

def SetTestSuite(suites):
    validation_suite = suites['validation']

    validation_suite.addTests(
         KratosUnittest.TestLoader().loadTestsFromTestCases([
          fluid_dem_coupling_one_way_test,
          sdem_pfem_coupling_one_way_test
          ])
    )

    return validation_suite

def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites
    night_suite = SetTestSuite(suites)
    suites['all'].addTests(night_suite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.DETAIL)
    KratosUnittest.runTests(AssembleTestSuites())
