from KratosMultiphysics import *
from KratosMultiphysics.RANSModellingApplication import *


def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    Tester.RunTestSuite("KratosRansFastSuite")


if __name__ == '__main__':
    run()
