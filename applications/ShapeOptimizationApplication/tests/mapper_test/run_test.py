# Import Kratos core and apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
from KratosMultiphysics.KratosUnittest import TestCase
from gid_output_process import GiDOutputProcess
import mapper_factory as mapper_factory
import math

# =======================================================================================================
# Auxiliary functions
# =======================================================================================================
def OutputResults(model_part, file_name):
    output_parameters = Parameters("""
    {
        "result_file_configuration" : {
            "gidpost_flags"       : {
                "GiDPostMode"           : "GiD_PostBinary",
                "WriteDeformedMeshFlag" : "WriteDeformed",
                "WriteConditionsFlag"   : "WriteConditions",
                "MultiFileFlag"         : "SingleFile"
            },
            "file_label"          : "step",
            "output_control_type" : "step",
            "output_frequency"    : 1,
            "nodal_results"       : ["CONTROL_POINT_UPDATE","CONTROL_POINT_CHANGE","SHAPE_UPDATE","PRESSURE","AIR_PRESSURE","WATER_PRESSURE"]
        },
        "point_data_configuration"  : []
    }""")
    gid_output_original = GiDOutputProcess(model_part, file_name, output_parameters)
    gid_output_original.ExecuteInitialize()
    gid_output_original.ExecuteBeforeSolutionLoop()
    gid_output_original.ExecuteInitializeSolutionStep()
    gid_output_original.PrintOutput()
    gid_output_original.ExecuteFinalizeSolutionStep()
    gid_output_original.ExecuteFinalize()

def Norm2OfVectorVariable(model_part, nodal_varible):
    norm_2 = 0
    for node in model_part.Nodes:
        temp_vec = node.GetSolutionStepValue(nodal_varible)
        norm_2 = norm_2 + temp_vec[0]*temp_vec[0] + temp_vec[1]*temp_vec[1] + temp_vec[2]*temp_vec[2]
    return math.sqrt(norm_2)

def Norm2OfScalarVariable(model_part, nodal_varible):
    norm_2 = 0
    for node in model_part.Nodes:
        temp_scalar = node.GetSolutionStepValue(nodal_varible)
        norm_2 = norm_2 + temp_scalar*temp_scalar
    return math.sqrt(norm_2)

# =======================================================================================================
# Set and read input data
# =======================================================================================================

# Import model parts
plate_with_trias = ModelPart("plate_with_trias")
plate_with_trias.AddNodalSolutionStepVariable(CONTROL_POINT_UPDATE)
plate_with_trias.AddNodalSolutionStepVariable(CONTROL_POINT_CHANGE)
plate_with_trias.AddNodalSolutionStepVariable(SHAPE_UPDATE)
plate_with_trias.AddNodalSolutionStepVariable(PRESSURE)
plate_with_trias.AddNodalSolutionStepVariable(WATER_PRESSURE)
plate_with_trias.AddNodalSolutionStepVariable(AIR_PRESSURE)
model_part_io = ModelPartIO("plate_with_trias")
model_part_io.ReadModelPart(plate_with_trias)

plate_with_quads = ModelPart("plate_with_quads")
plate_with_quads.AddNodalSolutionStepVariable(CONTROL_POINT_UPDATE)
plate_with_quads.AddNodalSolutionStepVariable(CONTROL_POINT_CHANGE)
plate_with_quads.AddNodalSolutionStepVariable(SHAPE_UPDATE)
model_part_io = ModelPartIO("plate_with_quads")
model_part_io.ReadModelPart(plate_with_quads)

# Set an input profile for the mapping variables (some saddle profile)
for node in plate_with_trias.Nodes:
    tmp_vector = [0.1, 0, (0.5-node.X)*(0.5-node.Y)]
    node.SetSolutionStepValue(CONTROL_POINT_UPDATE,tmp_vector)

for node in plate_with_trias.Nodes:
    node.SetSolutionStepValue(PRESSURE,(0.5-node.X)*(0.5-node.Y))

# =======================================================================================================
# Perform tests
# =======================================================================================================

# Test matrix-free mapper
mapper_settings = Parameters("""
{
    "filter_function_type"       : "linear",
    "filter_radius"              : 0.4,
    "max_nodes_in_filter_radius" : 10000,
    "matrix_free_filtering"      : true
}""")
matrix_mapper = mapper_factory.CreateMapper(plate_with_trias,plate_with_trias,mapper_settings)
matrix_mapper.Map(CONTROL_POINT_UPDATE,CONTROL_POINT_CHANGE)
matrix_mapper.InverseMap(CONTROL_POINT_CHANGE,SHAPE_UPDATE)

norm_2_result = Norm2OfVectorVariable(plate_with_trias, SHAPE_UPDATE)
TestCase().assertAlmostEqual(norm_2_result, 1.283132791556226, 12)

# Test matrix mapper
mapper_settings = Parameters("""
{
    "filter_function_type"       : "linear",
    "filter_radius"              : 0.4,
    "max_nodes_in_filter_radius" : 1000
}""")
matrix_mapper = mapper_factory.CreateMapper(plate_with_trias,plate_with_trias,mapper_settings)
matrix_mapper.Map(CONTROL_POINT_UPDATE,CONTROL_POINT_CHANGE)
matrix_mapper.InverseMap(CONTROL_POINT_CHANGE,SHAPE_UPDATE)

norm_2_result = Norm2OfVectorVariable(plate_with_trias, SHAPE_UPDATE)
TestCase().assertAlmostEqual(norm_2_result, 1.2831327915562258, 12)

# Test matrix mapper with consistent mapping
mapper_settings = Parameters("""
{
    "filter_function_type"       : "linear",
    "filter_radius"              : 0.4,
    "max_nodes_in_filter_radius" : 1000,
    "consistent_mapping"         : true
}""")
matrix_mapper = mapper_factory.CreateMapper(plate_with_trias,plate_with_trias,mapper_settings)
matrix_mapper.Map(CONTROL_POINT_UPDATE,CONTROL_POINT_CHANGE)
matrix_mapper.InverseMap(CONTROL_POINT_CHANGE,SHAPE_UPDATE)

norm_2_result = Norm2OfVectorVariable(plate_with_trias, SHAPE_UPDATE)
TestCase().assertAlmostEqual(norm_2_result, 1.266374348187224, 12)

# Test rectangular matrix mapper
mapper_settings = Parameters("""
{
    "filter_function_type"       : "linear",
    "filter_radius"              : 0.4,
    "max_nodes_in_filter_radius" : 1000
}""")
matrix_mapper = mapper_factory.CreateMapper(plate_with_trias,plate_with_quads,mapper_settings)
matrix_mapper.Map(CONTROL_POINT_UPDATE,CONTROL_POINT_CHANGE)
matrix_mapper.InverseMap(CONTROL_POINT_CHANGE,SHAPE_UPDATE)

norm_2_results_quad = Norm2OfVectorVariable(plate_with_quads, CONTROL_POINT_CHANGE)
norm_2_results_tria = Norm2OfVectorVariable(plate_with_trias, SHAPE_UPDATE)

TestCase().assertAlmostEqual(norm_2_results_quad, 2.5408880662655733, 12)
TestCase().assertAlmostEqual(norm_2_results_tria, 4.48736454850266, 12)

# Test rectangular matrix mapper with matrix free mapper
mapper_settings = Parameters("""
{
    "filter_function_type"       : "linear",
    "filter_radius"              : 0.4,
    "max_nodes_in_filter_radius" : 1000,
    "matrix_free_filtering"      : true
}""")
matrix_mapper = mapper_factory.CreateMapper(plate_with_trias,plate_with_quads,mapper_settings)
matrix_mapper.Map(CONTROL_POINT_UPDATE,CONTROL_POINT_CHANGE)
matrix_mapper.InverseMap(CONTROL_POINT_CHANGE,SHAPE_UPDATE)

norm_2_results_quad = Norm2OfVectorVariable(plate_with_quads, CONTROL_POINT_CHANGE)
norm_2_results_tria = Norm2OfVectorVariable(plate_with_trias, SHAPE_UPDATE)

TestCase().assertAlmostEqual(norm_2_results_quad, 2.5408880662655733, 12)
TestCase().assertAlmostEqual(norm_2_results_tria, 4.48736454850266, 12)

# Test improved integration
mapper_settings = Parameters("""
{
    "filter_function_type"       : "linear",
    "filter_radius"              : 0.4,
    "max_nodes_in_filter_radius" : 1000,
    "improved_integration"       : true,
    "integration_method"         : "gauss_integration",
    "number_of_gauss_points"     : 5
}""")
matrix_mapper = mapper_factory.CreateMapper(plate_with_trias,plate_with_trias,mapper_settings)
matrix_mapper.Map(CONTROL_POINT_UPDATE,CONTROL_POINT_CHANGE)
matrix_mapper.InverseMap(CONTROL_POINT_CHANGE,SHAPE_UPDATE)

norm_2_result = Norm2OfVectorVariable(plate_with_trias, SHAPE_UPDATE)
TestCase().assertAlmostEqual(norm_2_result, 1.3164625011428233, 12)

# Test scalar mapping
mapper_settings = Parameters("""
{
    "filter_function_type"       : "linear",
    "filter_radius"              : 0.4,
    "max_nodes_in_filter_radius" : 1000,
    "matrix_free_filtering"      : true
}""")
matrix_mapper = mapper_factory.CreateMapper(plate_with_trias,plate_with_trias,mapper_settings)
matrix_mapper.Map(PRESSURE,AIR_PRESSURE)
matrix_mapper.InverseMap(AIR_PRESSURE,WATER_PRESSURE)

norm_2_result = Norm2OfScalarVariable(plate_with_trias, WATER_PRESSURE)
TestCase().assertAlmostEqual(norm_2_result, 0.610521887077, 12)

# OutputResults(plate_with_trias,"results_tria_plate")
# OutputResults(plate_with_quads,"results_quad_plate")

# =======================================================================================================
# Clean folder
# =======================================================================================================

kratos_utilities.DeleteFileIfExisting("plate_with_trias.time")
kratos_utilities.DeleteFileIfExisting("plate_with_quads.time")

# =======================================================================================================
