from KratosMultiphysics import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)  # TESTS_OUTPUTS
    Tester.RunTestCases("*ResidualBasedAdjointBossak_TwoMassSpringDamperSystem_Conditions*")


if __name__ == '__main__':
    run()