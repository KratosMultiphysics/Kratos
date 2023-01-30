# Import Kratos core and apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.ShapeOptimizationApplication.analyzers.analyzer_base import AnalyzerBaseClass
from KratosMultiphysics.KratosUnittest import TestCase
import csv, os

# Read parameters
with open("parameters.json",'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

model = KM.Model()

class CustomAnalyzer(AnalyzerBaseClass):
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):
        mp = current_design.GetSubModelPart("edge_sens")
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


# Create optimizer and perform optimization
optimizer = optimizer_factory.Create(model, parameters["optimization_settings"], CustomAnalyzer())
optimizer.Optimize()

# =======================================================================================================
# Test results and clean directory
# =======================================================================================================
output_directory = parameters["optimization_settings"]["output"]["output_directory"].GetString()
optimization_log_filename = parameters["optimization_settings"]["output"]["optimization_log_filename"].GetString() + ".csv"
optimization_model_part_name = parameters["optimization_settings"]["model_settings"]["model_part_name"].GetString()

# Testing
with open(os.path.join(output_directory, optimization_log_filename), 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    last_line = None
    for line in reader:
        if not line:
            continue
        else:
            last_line = line

    resulting_optimization_iterations = int(last_line[0].strip())
    resulting_improvement = float(last_line[2].strip())

    # Check against specifications
    TestCase().assertEqual(resulting_optimization_iterations, 2)
    TestCase().assertAlmostEqual(resulting_improvement, -49.5142, 4)

expected_shape_change = {
    1: [-14.489666390559485, -8.365612827614228, 0.0],
    2: [10.48927962447755, -0.07307493673920817, 0.0],
    3: [-14.610503945724389, -8.435378422921769, 0.0],
    4: [10.711252869192506, -0.06583481797959133, 0.0],
    5: [-13.651092020770516, -7.66580088452228, 0.0],
    6: [-12.329971710549117, -6.39279549106124, 0.0],
    7: [-10.282755587645758, -4.868987133654232, 0.0],
    8: [-7.304385228960026, -3.489537463069345, 0.0],
    9: [-3.386710208733131, -2.3472493760049495, 0.0],
    10: [0.6737872321748706, -1.4479049341917931, 0.0],
    11: [4.094429978248886, -0.8243491932316478, 0.0],
    12: [6.8684518804921675, -0.42160320782139027, 0.0],
    13: [8.995917737587739, -0.1867540618978712, 0.0],
    14: [-13.698438205861017, -7.6920415965687745, 0.0],
    15: [-12.278779403746556, -6.3677441251508915, 0.0],
    16: [-10.124232426946033, -4.809772694184087, 0.0],
    17: [-7.042206147153214, -3.418824418997083, 0.0],
    18: [-3.0635964516130345, -2.282088405484788, 0.0],
    19: [0.990420969812569, -1.3937921370669817, 0.0],
    20: [4.3936799427993, -0.7831084092480116, 0.0],
    21: [7.145640644091602, -0.39324108908870103, 0.0],
    22: [9.247879015370504, -0.1702488119222248, 0.0],
}

_t = TestCase()
for node in model.GetModelPart(optimization_model_part_name).Nodes:
    shape_change = node.GetSolutionStepValue(KSO.SHAPE_CHANGE)
    _t.assertVectorAlmostEqual(shape_change, expected_shape_change[node.Id], 7)
    # print(f"{node.Id}: {[shape_change[0], shape_change[1], shape_change[2]]},")

# =======================================================================================================