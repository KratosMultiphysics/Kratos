from KratosMultiphysics import *

def runTests():
    Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
    Tester.RunTestCases("KratosCore*") # Is there a way to detect the cpp tests in the core?

if __name__ == '__main__':
    runTests()

