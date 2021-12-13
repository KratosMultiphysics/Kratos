import KratosMultiphysics
import test_WindEngineeringApplication

if __name__ == "__main__":
    if not KratosMultiphysics.IsDistributedRun():
        raise Exception("This test script can only be executed in MPI!")
    test_WindEngineeringApplication.Run(enable_mpi=True)