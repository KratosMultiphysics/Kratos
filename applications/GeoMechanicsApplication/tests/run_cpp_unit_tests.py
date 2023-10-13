import sys
import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication
import KratosMultiphysics.LinearSolversApplication

import argparse

def run():
    KratosMultiphysics.Tester.SetVerbosity(KratosMultiphysics.Tester.Verbosity.TESTS_LIST)

    exitCode1 = KratosMultiphysics.Tester.RunTestSuite("KratosGeoMechanicsFastSuite")

    parser = argparse.ArgumentParser()
    parser.add_argument("--onlyFastSuite", action='store_true')
    args = parser.parse_args()

    if args.onlyFastSuite :
        sys.exit(exitCode1)

    exitCode2 = KratosMultiphysics.Tester.RunTestSuite("KratosGeoMechanicsIntegrationSuite")

    if exitCode1 !=0 or exitCode2 !=0:
        print('\nAt least one test failed in one of the suites, please check your results!')

    sys.exit(exitCode1 if exitCode1 != 0 else exitCode2)

if __name__ == '__main__':
    run()
