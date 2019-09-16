"""
# IT'S ONLY A TEST. SCRIPT NOT COMPLETE 
"""


import KratosMultiphysics as Kratos
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.FluidDynamicsApplication

from geo_importer import GeoImporter
from geo_mesher import GeoMesher
from geo_preprocessor import GeoPreprocessor
from geo_building import GeoBuilding


# we read the JSON file for the CFD analysis
input_string = open("data/parameters/test/ProjectParameters.json",'r').read()
settings = KratosMultiphysics.Parameters(input_string)
print(settings.PrettyPrintJsonString()); input()
print("01")

# we set "problem_name"
problem_name = "test_00"
settings["problem_data"].AddEmptyValue("problem_name")
settings["problem_data"]["problem_name"].SetString(problem_name)

print("02")

# we set "output_name"
problem_name = "ABC"
settings["output_processes"]["gid_output"]["Parameters"].AddEmptyValue("output_name")
settings["output_processes"]["gid_output"]["Parameters"]["output_name"].SetString(problem_name)

print("03")

# we set "model_import_settings"
input_filename = "10_Box_buildings_subtracted_3_after_CleanConditions"
settings["solver_settings"]["model_import_settings"].AddEmptyValue("input_filename")
settings["solver_settings"]["model_import_settings"]["input_filename"].SetString(input_filename)

print("04")

# we set "volume_model_part_name"
volume_name = "Parts_Fluid"
settings["solver_settings"].AddEmptyValue("volume_model_part_name")
settings["solver_settings"]["volume_model_part_name"].SetString(volume_name)
# # we set "skin_parts"
# skin_name = ["Inlet", "Outlet", "BottomModelPart", "TopModelPart", "SKIN_ISOSURFACE"]
# settings["solver_settings"].AddEmptyValue("skin_parts")
# settings["solver_settings"]["skin_parts"].SetString(skin_name)

print("05")

# we set sub mdoel part name
sub_model_inlet = "FluidModelPart.Inlet"
settings["processes"]["boundary_conditions_process_list"]["Parameters"].AddEmptyValue("model_part_name")
settings["processes"]["boundary_conditions_process_list"]["Parameters"]["model_part_name"].SetString(sub_model_inlet)

