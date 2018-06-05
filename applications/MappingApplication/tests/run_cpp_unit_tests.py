from KratosMultiphysics import *
try: # Importing MPI before mapping makes the logo only print once
    import KratosMultiphysics.mpi as KratosMPI
except:
    pass
from KratosMultiphysics.MappingApplication import *

import KratosMultiphysics.kratos_utilities as kratos_utils

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
            Logger.PrintInfo("\ncpp tests MappingApplication", "Creating the partitioned files\n")

        import create_partitioned_mdpa_utility
        num_partitions = KratosMPI.mpi.size
        create_partitioned_mdpa_utility.CreatePartitionedMdpaFiles(num_partitions)
        # Tester.RunTestSuite("KratosMappingApplicationMPITestSuite")

    # This suite contains tests that work both with and without MPI
    # Note that if the partitioned mdpa files were not created before the MPI tests are skipped
    Logger.PrintInfo("\ncpp tests MappingApplication", "Running General tests\n")
    Tester.RunTestSuite("KratosMappingApplicationGeneralTestSuite")

    folder_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "partitioned_mdpa_files")
    kratos_utils.DeleteDirectoryIfExisting(folder_path)

if __name__ == '__main__':
    run()
