from KratosMultiphysics import *
from KratosMultiphysics.StatisticsApplication import *


def run():
    Tester.SetVerbosity(Tester.Verbosity.FAILED_TESTS_OUTPUTS)  # TESTS_OUTPUTS
    Tester.RunTestSuite("KratosStatisticsFastSuite")


if __name__ == '__main__':
    run()