from KratosMultiphysics import *
from KratosMultiphysics.MappingApplication import *

def IsMPIExecution():
    try:
        import KratosMultiphysics.mpi as KratosMPI
        print("here", KratosMPI.mpi.size)
        if KratosMPI.mpi.size > 1:
            return True
        else:
            return False
    except:
        return False

def run():
    Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
    if not IsMPIExecution():
        # This suite contains tests that either don't run with MPI
        # or where it doesn't make sense to run in MPI
        Tester.RunTestSuite("KratosMappingApplicationSerialTestSuite")
    else:
        # This suite contains tests that require MPI
        Tester.RunTestSuite("KratosMappingApplicationMPITestSuite")

    # This suite runs both with and without MPI
    Tester.RunTestSuite("KratosMappingApplicationGeneralTestSuite")

if __name__ == '__main__':
    run()
