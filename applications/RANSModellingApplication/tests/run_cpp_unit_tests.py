from KratosMultiphysics import *
from KratosMultiphysics.RANSModellingApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.FAILED_TESTS_OUTPUTS)
    Tester.RunTestSuite("RansEvmKEpsilonConditionInterfaces")
    Tester.RunTestSuite("RansEvmKEpsilonConditionMethods")
    Tester.RunTestSuite("RansEvmKEpsilonElementInterfaces")
    Tester.RunTestSuite("RansEvmKEpsilonElementMethods")

if __name__ == '__main__':
    run()
