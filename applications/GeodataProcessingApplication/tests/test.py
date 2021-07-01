import os
import time
import json

# Kratos import
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.FluidDynamicsApplication as KratosFluidDynamics
import KratosMultiphysics.MeshingApplication as KratosMeshing
import KratosMultiphysics.MappingApplication as KratosMapping


# Kratos Fluid Dynamic Analysis Imports
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
import KratosMultiphysics.TrilinosApplication.trilinos_linear_solver_factory as trilinos_linear_solver_factory

from KratosMultiphysics.mpi.distributed_gid_output_process import DistributedGiDOutputProcess as GiDOutputProcess
import KratosMultiphysics.mpi.distributed_import_model_part_utility as distributed_import_model_part_utility


# GeodataProcessing import
from KratosMultiphysics.GeodataProcessingApplication.geo_processor import GeoProcessor
from KratosMultiphysics.GeodataProcessingApplication.geo_importer import GeoImporter
from KratosMultiphysics.GeodataProcessingApplication.geo_mesher import GeoMesher
from KratosMultiphysics.GeodataProcessingApplication.geo_preprocessor import GeoPreprocessor
from KratosMultiphysics.GeodataProcessingApplication.geo_building import GeoBuilding
from KratosMultiphysics.GeodataProcessingApplication.geo_model import GeoModel
from KratosMultiphysics.GeodataProcessingApplication.geo_data import GeoData



def _OutputToGid(main_model_part, file_name):
    gid_output = GiDOutputProcess(
        main_model_part,
        "gid_output/"+file_name,
        KratosMultiphysics.Parameters("""
            {
                "result_file_configuration" : {
                    "gidpost_flags": {
                        "GiDPostMode": "GiD_PostBinary",
                        "WriteConditionsFlag": "WriteConditions",
                        "MultiFileFlag": "SingleFile"
                    },
                    "nodal_results"       : ["DISTANCE", "PARTITION_INDEX"],
                    "nodal_nonhistorical_results": [],
                    "gauss_point_results" : [],
                    "nodal_flags_results": [],
                    "elemental_conditional_flags_results": []
                }
            }
            """)
    )
    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()



with open("data/json/ProjectParameters.json",'r') as parameter_file:
    parameters = KratosMultiphysics.Parameters(parameter_file.read())
    print("parameters read!")

with open("data/json/RemeshingParameters.json",'r') as parameter_file:
    remeshing_parameters = KratosMultiphysics.Parameters(parameter_file.read())
    print("remeshing_parameters read!")


importer = GeoImporter()

# import domain from mdpa file
importer._InitializeModelPart("MainModelPart")
importer.ModelPart.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
main_model_part = importer.ModelPart
model_part_in = parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()
# KratosMultiphysics.ModelPartIO(model_part_in).ReadModelPart(main_model_part)


import_parameters = KratosMultiphysics.Parameters("""{
    "echo_level" : 0,
    "model_import_settings" : {
        "input_type" : "mdpa",
        "input_filename" : "model_part_name",
        "partition_in_memory" : false
    }
}""")

import_parameters["model_import_settings"]["input_filename"].SetString(parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString())



# Construct the Trilinos import model part utility
distributed_model_part_importer = KratosMultiphysics.mpi.distributed_import_model_part_utility.DistributedImportModelPartUtility(main_model_part, import_parameters)


# Execute the Metis partitioning and reading
distributed_model_part_importer.ImportModelPart()
distributed_model_part_importer.CreateCommunicators()

communicator = KratosMultiphysics.DataCommunicator.GetDefault()
rank = communicator.Rank()
size = communicator.Size()

# if rank == 0:
#     print(import_parameters)
#     input()
# #     print(main_model_part)
# communicator.Barrier()

if rank == 0:
    if not os.path.exists("gid_output"):
        os.makedirs("gid_output")

communicator.Barrier()


local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(
                main_model_part,
                KratosMultiphysics.DISTANCE,
                KratosMultiphysics.DISTANCE_GRADIENT,
                KratosMultiphysics.NODAL_AREA)
local_gradient.Execute()

find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
find_nodal_h.Execute()

metric_parameters = remeshing_parameters["metric_parameters"]
metric_process = KratosMultiphysics.MeshingApplication.ComputeLevelSetSolMetricProcess3D(
                main_model_part,
                KratosMultiphysics.DISTANCE_GRADIENT,
                metric_parameters)
metric_process.Execute()

parmmg_parameters = remeshing_parameters["parmmg_parameters"]
parmmg_process = KratosMultiphysics.MeshingApplication.ParMmgProcess3D(
                main_model_part,
                parmmg_parameters)
parmmg_process.Execute()

_OutputToGid(main_model_part, 'final_remeshed_distance_field')

