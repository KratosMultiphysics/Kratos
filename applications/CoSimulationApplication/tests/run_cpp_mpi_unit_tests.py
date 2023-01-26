import KratosMultiphysics as KM

if not KM.IsDistributedRun():
    raise Exception("This test script can only be executed in MPI!")

from KratosMultiphysics.CoSimulationApplication import * # registering tests

def run():
    KM.Tester.SetVerbosity(KM.Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    KM.Tester.RunTestSuite("KratosCosimulationMPIFastSuite")

if __name__ == '__main__':
    run()
