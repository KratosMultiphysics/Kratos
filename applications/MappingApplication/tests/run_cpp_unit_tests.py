from KratosMultiphysics import *
try: # Importing MPI before mapping maked the logo only print once
    import KratosMultiphysics.mpi as KratosMPI
except:
    pass
from KratosMultiphysics.MappingApplication import *

def IsMPIExecution():
    try:
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
        Logger.PrintInfo("\ncpp tests MappingApplication", "Running Serial tests\n")
        Tester.RunTestSuite("KratosMappingApplicationSerialTestSuite")
    else:
        # This suite contains tests that require MPI
        if KratosMPI.mpi.rank == 0:
            Logger.PrintInfo("\ncpp tests MappingApplication", "Running MPI tests\n")
        Tester.RunTestSuite("KratosMappingApplicationMPITestSuite")

    # This suite runs both with and without MPI
    # Tester.RunTestSuite("KratosMappingApplicationGeneralTestSuite")

if __name__ == '__main__':
    run()
