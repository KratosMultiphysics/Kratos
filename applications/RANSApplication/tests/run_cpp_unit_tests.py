from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.RANSApplication import *


def run():
    print(dir(Tester.Verbosity))
    Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)  # TESTS_OUTPUTS
    Tester.RunTestCases("*TestRansKEpsilonEpsilonCWD2D3N_CalculateLocalVelocityCon*")



if __name__ == '__main__':
    run()
