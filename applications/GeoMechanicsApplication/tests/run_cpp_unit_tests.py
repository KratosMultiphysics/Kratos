import KratosMultiphysics
from KratosMultiphysics.GeoMechanicsApplication import *


def run():
    KratosMultiphysics.Tester.SetVerbosity(KratosMultiphysics.Tester.Verbosity.TESTS_LIST)
    return KratosMultiphysics.Tester.RunTestSuite("KratosGeoMechanicsFastSuite")


if __name__ == '__main__':
    run()
