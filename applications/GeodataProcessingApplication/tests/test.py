# test mpi

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.MeshingApplication as KratosMesh


if (Kratos.IsDistributedRun()):
    import KratosMultiphysics.mpi as KratosMPI
    print("[DEBUG] Kratos Parallel")
    
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()


# import sys
# sys.path.append('../python_scripts/')

# from GeodataProcessingApplication.python_scripts.geo_importer import GeoImporter
from python_scripts.geo_importer import GeoImporter

# from geo_importer import GeoImporter
# from geo_mesher import GeoMesher
# from geo_preprocessor import GeoPreprocessor
# from geo_building import GeoBuilding
# from geo_model import GeoModel
# from geo_data import GeoData

import math
import os
import time

