# Definition of the classes for the NIGHTLY TESTS

#Iimport Kratos
import KratosMultiphysics
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.SwimmingDEMApplication
from KratosMultiphysics import Logger

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import SmallTests
import GentleInjectionAndErasureTestFactory as GentleTF
import PeriodicBackwardCouplingTestFactory as PeriodicBackwardTF


# class gentle_injection_test(GentleTF.GentleInjectionAndErasureTestFactory):
#      file_name = "fluid_convergence_tests/cube_cavity_with_inlet"
#      file_parameters_harsh = "fluid_convergence_tests/ProjectParametersInjectionHarsh.json"
#      file_parameters_gentle = "fluid_convergence_tests/ProjectParametersInjectionGentle.json"
#      def GetGentleParameterValueAndName(self, parameters):
#          parameter_name = 'initiation_interval'
#          return parameters['coupling']['gentle_coupling_initiation'][parameter_name].GetDouble(), parameter_name
# class gentle_erasure_test(GentleTF.GentleInjectionAndErasureTestFactory):
#      file_name = "fluid_convergence_tests/cube_cavity"
#      file_parameters_harsh = "fluid_convergence_tests/ProjectParametersErasureHarsh.json"
#      file_parameters_gentle = "fluid_convergence_tests/ProjectParametersErasureGentle.json"
#      def GetGentleParameterValueAndName(self, parameters):
#          parameter_name = 'destruction_delay_interval'
#          return parameters['dem_parameters']['creator_destructor_settings'][parameter_name].GetDouble(), parameter_name

# class periodic_backward_coupling(PeriodicBackwardTF.PeriodicBackwardCouplingTestFactory):
#      file_name = "fluid_convergence_tests/cube_cavity"
#      file_parameters_harsh = "fluid_convergence_tests/ProjectParametersPeriodicHarsh.json"
#      file_parameters_gentle = "fluid_convergence_tests/ProjectParametersPeriodicHarsh.json"
#      def GetGentleParameterValueAndName(self, parameters):
#          parameter_name = 'initiation_interval'
#          return parameters['coupling']['gentle_coupling_initiation'][parameter_name].GetDouble(), parameter_name

# class gentle_injection_test(GentleTF.GentleInjectionAndErasureTestFactory):
#      file_name = "fluid_convergence_tests/cube_cavity_with_inlet"
#      file_parameters_harsh = "fluid_convergence_tests/ProjectParametersInjectionHarshPeriodic.json"
#      file_parameters_gentle = "fluid_convergence_tests/ProjectParametersInjectionHarshPeriodic.json"
#      def GetGentleParameterValueAndName(self, parameters):
#          parameter_name = 'initiation_interval'
#          return parameters['coupling']['gentle_coupling_initiation'][parameter_name].GetDouble

# List of tests that are available
available_tests = []
available_tests += [test_class for test_class in GentleTF.GentleInjectionAndErasureTestFactory.__subclasses__()]
available_tests += [test_class for test_class in PeriodicBackwardTF.PeriodicBackwardCouplingTestFactory.__subclasses__()]

def SetTestSuite(suites):
    night_suite = suites['nightly']

    night_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(available_tests))

    return night_suite

def AssembleTestSuites():
    suites = KratosUnittest.KratosSuites
    # small_suite = SmallTests.SetTestSuite(suites)
    # suites['all'].addTests(small_suite)
    night_suite = SetTestSuite(suites)
    suites['all'].addTests(night_suite)

    return suites

if __name__ == '__main__':
    debug_mode = False
    if debug_mode:
        severity = Logger.Severity.DETAIL
    else:
        severity = Logger.Severity.WARNING
    Logger.GetDefaultOutput().SetSeverity(severity)
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())
