from KratosMultiphysics import *
from KratosMultiphysics.RANSModellingApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.FAILED_TESTS_OUTPUTS) # TESTS_OUTPUTS
    Tester.RunTestSuite("RANSYPlusModelSensitivities")
    Tester.RunTestSuite("RANSModellingApplicationElementInterfaces")
    Tester.RunTestSuite("RANSEvModelsKEpsilonElementMethods")
    Tester.RunTestSuite("RANSEvModelsKEpsilonConditionMethods")
    Tester.RunTestSuite("RANSModellingApplicationConditionInterfaces")

if __name__ == '__main__':
    run()
