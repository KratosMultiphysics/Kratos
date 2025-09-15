# Import Kratos core and apps
import KratosMultiphysics as Kratos

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.ShapeOptimizationApplication.analyzers.analyzer_base import AnalyzerBaseClass
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess

# Auxiliary imports
# from KratosMultiphysics.vtk_output_process import VtkOutputProcess
# from KratosMultiphysics.json_output_process import JsonOutputProcess

# =======================================================================================================
# Auxiliary functions
# =======================================================================================================
# def OutputResults(model):
#     output_parameters = Kratos.Parameters("""
#      {
#         "model_part_name"        : "hyperbolic_paraboloid.design_surface",
#         "write_ids"              : false,
#         "file_format"            : "binary",
#         "output_sub_model_parts" : false,
#         "output_path"            : "Optimization_Results",
#         "nodal_solution_step_data_variables": ["MAX_NEIGHBOUR_DISTANCE", "GAUSSIAN_CURVATURE","VERTEX_MORPHING_RADIUS_RAW","VERTEX_MORPHING_RADIUS", "DF1DX","DF1DX_MAPPED", "CONTROL_POINT_UPDATE", "SHAPE_UPDATE"]
#     }
#     """)
#     vtk_output_original = VtkOutputProcess(model, output_parameters)
#     vtk_output_original.ExecuteInitialize()
#     vtk_output_original.ExecuteBeforeSolutionLoop()
#     vtk_output_original.ExecuteInitializeSolutionStep()
#     vtk_output_original.PrintOutput()
#     vtk_output_original.ExecuteFinalizeSolutionStep()
#     vtk_output_original.ExecuteFinalize()

# =======================================================================================================
# Set and read input data
# =======================================================================================================

# Read parameters
with open("optimization_parameters.json",'r') as parameter_file:
    parameters = Kratos.Parameters(parameter_file.read())


class CustomAnalyzer(AnalyzerBaseClass):
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):
        if communicator.isRequestingValueOf("no_objective"):
            value = 0.0
            communicator.reportValue("no_objective", value)

        if communicator.isRequestingGradientOf("no_objective"):
            gradient = {}
            for node in current_design.Nodes:
                if node.X == -1 or node.X == 1:
                    gradient[node.Id] = [0.0, 0.0, -1.0]
                if node.Y == -1 or node.Y == 1:
                    if node.X == -1 or node.X == 1:
                        gradient[node.Id] = [0.0, 0.0, 0.0]
                    else:
                        gradient[node.Id] = [0.0, 0.0, 1.0]
            communicator.reportGradient("no_objective", gradient)

model = Kratos.Model()

# Create optimizer and perform optimization
optimizer = optimizer_factory.Create(model, parameters["optimization_settings"], CustomAnalyzer())
optimizer.Optimize()

optimization_model_part_name = parameters["optimization_settings"]["model_settings"]["model_part_name"].GetString()

# =======================================================================================================
# Perform tests
# =======================================================================================================
# # write json output
# output_process = JsonOutputProcess(model, Kratos.Parameters(
#     """{
#         "output_variables" : ["DF1DX_MAPPED", "SHAPE_UPDATE"],
#         "output_file_name" : "gaussian_curvature_results.json",
#         "model_part_name"  : \""""+optimization_model_part_name+"""\",
#         "time_frequency"   : 0.0
#     }"""))

# output_process.ExecuteInitialize()
# output_process.ExecuteBeforeSolutionLoop()
# output_process.ExecuteInitializeSolutionStep()
# output_process.ExecuteFinalizeSolutionStep()
# output_process.ExecuteFinalize()

check_process = FromJsonCheckResultProcess(model, Kratos.Parameters(
    """{
        "check_variables"  : ["DF1DX_MAPPED", "SHAPE_UPDATE"],
        "input_file_name"  : "gaussian_curvature_results.json",
        "model_part_name"  : \""""+optimization_model_part_name+"""\",
        "time_frequency"   : 0.0
    }"""))
check_process.ExecuteInitialize()
check_process.ExecuteBeforeSolutionLoop()
check_process.ExecuteInitializeSolutionStep()
check_process.ExecuteFinalizeSolutionStep()
check_process.ExecuteFinalize()
# =======================================================================================================
