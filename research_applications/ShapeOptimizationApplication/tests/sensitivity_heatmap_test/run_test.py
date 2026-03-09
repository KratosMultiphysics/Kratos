# Import Kratos core and apps
import KratosMultiphysics as Kratos
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.ShapeOptimizationApplication.analyzers.analyzer_base import AnalyzerBaseClass
from KratosMultiphysics.ShapeOptimizationApplication.utilities import custom_math as cm
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess
# from KratosMultiphysics.json_output_process import JsonOutputProcess

# Read parameters
with open("parameters.json",'r') as parameter_file:
    parameters = Kratos.Parameters(parameter_file.read())

model = Kratos.Model()

class CustomAnalyzer(AnalyzerBaseClass):
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):
        mp = current_design.GetSubModelPart("design_surface")
        if communicator.isRequestingValueOf("x_squared_sum"):
            value = 0.0
            for node in mp.Nodes:
                value += node.X**2
            communicator.reportValue("x_squared_sum", value)

        if communicator.isRequestingGradientOf("x_squared_sum"):
            gradient = {node.Id: [0.0, 0.0, 0.0] for node in current_design.Nodes}
            for node in mp.Nodes:
                gradient[node.Id] = [2*node.X, 0.0, 0.0]
            communicator.reportGradient("x_squared_sum", gradient)

        if communicator.isRequestingValueOf("shape_change"):
            value = 0.0
            for node in mp.Nodes:
                shape_change = node.GetSolutionStepValue(KSO.SHAPE_CHANGE)
                norm = cm.Norm2(shape_change)
                if norm > 2.0:
                    value += (norm - 2)**2
            communicator.reportValue("shape_change", value)

        if communicator.isRequestingGradientOf("x_squared_sum"):
            gradient = {node.Id: [0.0, 0.0, 0.0] for node in current_design.Nodes}
            for node in mp.Nodes:
                shape_change = node.GetSolutionStepValue(KSO.SHAPE_CHANGE)
                norm = cm.Norm2(shape_change)
                if norm > 2.0:
                    gradient[node.Id] = 2*(norm - 2) * shape_change
            communicator.reportGradient("shape_change", gradient)


# Create optimizer and perform optimization
optimizer = optimizer_factory.Create(model, parameters["optimization_settings"], CustomAnalyzer())
optimizer.Optimize()

# =======================================================================================================
# Test results and clean directory
# =======================================================================================================
output_directory = parameters["optimization_settings"]["output"]["output_directory"].GetString()
optimization_log_filename = parameters["optimization_settings"]["output"]["optimization_log_filename"].GetString() + ".csv"
optimization_model_part_name = parameters["optimization_settings"]["model_settings"]["model_part_name"].GetString()

# write json output
# output_process = JsonOutputProcess(model, Kratos.Parameters(
#     """{
#         "output_variables" : ["HEATMAP_DF1DX_X","HEATMAP_DF1DX_Y","HEATMAP_DF1DX_Z",
#                               "HEATMAP_DC1DX_X","HEATMAP_DC1DX_Y","HEATMAP_DC1DX_Z",
#                               "HEATMAP_L2"],
#         "output_file_name" : "heatmap_results.json",
#         "model_part_name"  : \""""+optimization_model_part_name+"""\",
#         "time_frequency"   : 0.0
#     }"""))

# output_process.ExecuteInitialize()
# output_process.ExecuteBeforeSolutionLoop()
# output_process.ExecuteInitializeSolutionStep()
# output_process.ExecuteFinalizeSolutionStep()
# output_process.ExecuteFinalize()

# Testing
check_process = FromJsonCheckResultProcess(model, Kratos.Parameters(
    """{
        "check_variables"  : ["HEATMAP_DF1DX_X","HEATMAP_DF1DX_Y","HEATMAP_DF1DX_Z",
                              "HEATMAP_DC1DX_X","HEATMAP_DC1DX_Y","HEATMAP_DC1DX_Z",
                              "HEATMAP_L2"],
        "input_file_name"  : "heatmap_results.json",
        "model_part_name"  : \""""+optimization_model_part_name+"""\",
        "time_frequency"   : 0.0
    }"""))
check_process.ExecuteInitialize()
check_process.ExecuteBeforeSolutionLoop()
check_process.ExecuteInitializeSolutionStep()
check_process.ExecuteFinalizeSolutionStep()
check_process.ExecuteFinalize()

# =======================================================================================================