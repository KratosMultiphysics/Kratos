from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
def import_solver(SolverSettings):
    """this function imports a solver named "solver_type" from SolverSettings
    solver_type is expected to be the FILENAME of the solver to be imported"""
    obj = __import__(SolverSettings.solver_type)
    return obj

def IsRankZero():
    try:
        import KratosMultiphysics.mpi as KratosMPI
        KratosMPI.mpi.world.barrier()
        my_rank = KratosMPI.mpi.rank
    except ImportError:
        my_rank = 0

    return my_rank == 0

def DeleteFileIfExisting(file_name):
    if IsRankZero():
        import os
        if os.path.isfile(file_name):
            os.remove(file_name)

def DeleteDirectoryIfExisting(directory_name):
    if IsRankZero():
        import os
        import shutil
        if os.path.isdir(directory_name):
            shutil.rmtree(directory_name)
