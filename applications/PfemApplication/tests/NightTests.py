# Definition of the classes for the NIGHTLY TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities

class newtonian_sloshing_2D(TF.TestFactory):
     file_name = "fluid_tests/newtonian/sloshing_2D"
     file_parameters = None
class newtonian_sloshing_3D(TF.TestFactory):
     file_name = "fluid_tests/newtonian/sloshing_3D"
     file_parameters = None

def SetTestSuite(suites):
     night_suite = suites['nightly']

     has_pfem_fluid_application = kratos_utilities.CheckIfApplicationsAvailable("PfemFluidMechanicsApplication")
 
     if ( has_pfem_fluid_application):
        night_suite.addTests(
             KratosUnittest.TestLoader().loadTestsFromTestCases([
                  newtonian_sloshing_2D
             ])
        )

     return night_suite
