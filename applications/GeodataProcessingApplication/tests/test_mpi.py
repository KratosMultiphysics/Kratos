# test mpi

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.MeshingApplication as KratosMesh

from KratosMultiphysics.testing.utilities import ReadModelPart

if (Kratos.IsDistributedRun()):
    import KratosMultiphysics.mpi as KratosMPI
    print("[DEBUG] Kratos Parallel")
    
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

from KratosMultiphysics.GeodataProcessingApplication.geo_processor import GeoProcessor
from KratosMultiphysics.GeodataProcessingApplication.geo_importer import GeoImporter
from KratosMultiphysics.GeodataProcessingApplication.geo_mesher import GeoMesher
from KratosMultiphysics.GeodataProcessingApplication.geo_preprocessor import GeoPreprocessor
from KratosMultiphysics.GeodataProcessingApplication.geo_building import GeoBuilding
from KratosMultiphysics.GeodataProcessingApplication.geo_model import GeoModel
from KratosMultiphysics.GeodataProcessingApplication.geo_data import GeoData

import math
import os
import time


def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName


num_test = "038_mpi"
print("\n\n[DEBUG] TEST ", num_test)
print("[DEBUG][start]", time.strftime("%H:%M:%S", time.localtime()))

if rank == 0:
    # we create a new folders for this test
    if not os.path.exists("cfd_data/test_{}".format(num_test)):
        os.mkdir("cfd_data/test_{}".format(num_test))
    if not os.path.exists("cfd_data/test_{}/gid_file".format(num_test)):
        os.mkdir("cfd_data/test_{}/gid_file".format(num_test))
    if not os.path.exists("cfd_data/test_{}/gid_file/Binary".format(num_test)):
        os.mkdir("cfd_data/test_{}/gid_file/Binary".format(num_test))
    if not os.path.exists("cfd_data/test_{}/gid_file/Ascii".format(num_test)):
        os.mkdir("cfd_data/test_{}/gid_file/Ascii".format(num_test))
    if not os.path.exists("cfd_data/test_{}/stl_file".format(num_test)):
        os.mkdir("cfd_data/test_{}/stl_file".format(num_test))
    if not os.path.exists("cfd_data/test_{}/analysis_file".format(num_test)):
        os.mkdir("cfd_data/test_{}/analysis_file".format(num_test))
    if not os.path.exists("cfd_data/test_{}/mdpa_file".format(num_test)):
        os.mkdir("cfd_data/test_{}/mdpa_file".format(num_test))

    if not os.path.exists("cfd_data/test_{}/AsterGDEM".format(num_test)):
        os.mkdir("cfd_data/test_{}/AsterGDEM".format(num_test))
    if not os.path.exists("cfd_data/test_{}/OSM".format(num_test)):
        os.mkdir("cfd_data/test_{}/OSM".format(num_test))

importer = GeoImporter()

# import domain from mdpa file
importer._InitializeModelPart("MainModelPart")
main_model_part = importer.ModelPart
main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(main_model_part)

# model_part_in = "data/mdpa_file/test_mpi/01_Mesh_cylinder"
model_part_in = "data/mdpa_file/test_mpi/02_Mesh_cylinder_distance_field_1"
KratosMultiphysics.ModelPartIO(model_part_in).ReadModelPart(main_model_part)

num_nodes = main_model_part.NumberOfNodes()

# for node in main_model_part.Nodes:
for i_n, node in enumerate(main_model_part.Nodes):
    if i_n < round(num_nodes/size + 10**-9):
        node.SetSolutionStepValue(Kratos.PARTITION_INDEX, 0)
    elif i_n < round(num_nodes/size * 2 + 10**-9):
        node.SetSolutionStepValue(Kratos.PARTITION_INDEX, 1)
    elif i_n < round(num_nodes/size * 3 + 10**-9):
        node.SetSolutionStepValue(Kratos.PARTITION_INDEX, 2)
    else:
        node.SetSolutionStepValue(Kratos.PARTITION_INDEX, 3)

if rank == 0:
    for i_n, node in enumerate(main_model_part.Nodes):
        print(node.GetSolutionStepValue(Kratos.PARTITION_INDEX))
"""
    DIVIDERE I NODI IN GRUPPI
"""


communicator = main_model_part.GetCommunicator().GetDataCommunicator()

local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(main_model_part,
                                                                  KratosMultiphysics.DISTANCE,
                                                                  KratosMultiphysics.DISTANCE_GRADIENT,
                                                                  KratosMultiphysics.NODAL_AREA)

local_gradient.Execute()

find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
find_nodal_h.Execute()

# We define a metric using the ComputeLevelSetSolMetricProcess
min_size = 10.0
max_size = 500.0
max_dist = 1.0
levelset_parameters = KratosMultiphysics.Parameters("""
    {
        "minimal_size"                         : """ + str(min_size) + """,
        "maximal_size"                         : """ + str(max_size) + """,
        "sizing_parameters":
        {
            "reference_variable_name"          : "DISTANCE",
            "boundary_layer_max_distance"      : """ + str(max_dist) + """,
            "interpolation"                    : "linear"
        },
        "enforce_current"                      : true,
        "anisotropy_remeshing"                 : true
    }
    """)

metric_process = KratosMesh.ComputeLevelSetSolMetricProcess3D(main_model_part,
                                                              KratosMultiphysics.DISTANCE_GRADIENT,
                                                              levelset_parameters)
metric_process.Execute()


pmmg_parameters = KratosMultiphysics.Parameters("""
    {
        "filename"                         : "output"
    }
    """)
pmmg_parameters["filename"].SetString(GetFilePath(pmmg_parameters["filename"].GetString()))

# print(main_model_part)
# input("vedi model part")

KratosMPI.ParallelFillCommunicator(main_model_part).Execute()

pmmg_process = KratosMultiphysics.MeshingApplication.ParMmgProcess3D(main_model_part.GetRootModelPart(), pmmg_parameters)
print("PARMMG")
pmmg_process.Execute()

# GiD file
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/02_Mesh_cylinder_distance_field_1".format(num_test), "GiD_PostAscii")
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/02_Mesh_cylinder_distance_field_1".format(num_test), "GiD_PostBinary")
# MDPA file
mesher.WriteMdpaOutput("cfd_data/test_{}/mdpa_file/02_Mesh_cylinder_distance_field_1".format(num_test))
