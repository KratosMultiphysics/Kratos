from __future__ import print_function, absolute_import, division

# Import Kratos
import KratosMultiphysics
try:
    import KratosMultiphysics.mpi as KratosMPI
    rank = KratosMPI.mpi.rank
except:
    rank = 0

# Other imports
import os

def DeleteFileIfExisting(file_name):
    if rank == 0:
        try:
            os.remove(file_name)
        except FileNotFoundError:
            pass
