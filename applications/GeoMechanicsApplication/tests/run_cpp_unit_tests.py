import sys
import KratosMultiphysics
from KratosMultiphysics.GeoMechanicsApplication import *


def run():
    KratosMultiphysics.Tester.SetVerbosity(KratosMultiphysics.Tester.Verbosity.TESTS_LIST)
    exitCode = KratosMultiphysics.Tester.RunTestSuite("KratosGeoMechanicsFastSuite")
    sys.exit(exitCode)


if __name__ == '__main__':
    run()
