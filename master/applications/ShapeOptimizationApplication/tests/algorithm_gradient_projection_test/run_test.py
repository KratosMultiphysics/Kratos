# Import Kratos core and apps
import KratosMultiphysics as KM

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication.analyzers.analyzer_base import AnalyzerBaseClass
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.KratosUnittest import TestCase
import KratosMultiphysics.kratos_utilities as kratos_utilities
import csv, os

# =======================================================================================================
# Define external analyzer
# =======================================================================================================
class Target():

    target_z = 2.0

    @classmethod
    def Value(cls, model_part):
        f = 0
        for node in model_part.Nodes:
            f += (cls.target_z - node.Z) ** 2
        return f

    @classmethod
    def Gradient(cls, model_part):
        gradient = {}
        for node in model_part.Nodes:
            sz = - 2 * (cls.target_z - node.Z)
            gradient[node.Id] = [0.0, 0.0, sz]
        return gradient

class CoordZ():

    @staticmethod
    def Value(model_part, node_id):
        return model_part.Nodes[node_id].Z

    @staticmethod
    def Gradient(model_part, node_id):
        gradient = {}
        for node in model_part.Nodes:
            if node.Id == node_id:
                gradient[node.Id] = [0.0, 0.0, 1.0]
            else:
                gradient[node.Id] = [0.0, 0.0, 0.0]
        return gradient


class CustomAnalyzer(AnalyzerBaseClass):

    # --------------------------------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):

        if communicator.isRequestingValueOf("target"):
            communicator.reportValue("target", Target.Value(current_design))
        if communicator.isRequestingGradientOf("target"):
            communicator.reportGradient("target", Target.Gradient(current_design))

        if communicator.isRequestingValueOf("z_15"):
            communicator.reportValue("z_15", CoordZ.Value(current_design, 15))
        if communicator.isRequestingGradientOf("z_15"):
            communicator.reportGradient("z_15", CoordZ.Gradient(current_design, 15))

        if communicator.isRequestingValueOf("z_40"):
            communicator.reportValue("z_40", CoordZ.Value(current_design, 40))
        if communicator.isRequestingGradientOf("z_40"):
            communicator.reportGradient("z_40", CoordZ.Gradient(current_design, 40))

        if communicator.isRequestingValueOf("z_65"):
            communicator.reportValue("z_65", CoordZ.Value(current_design, 65))
        if communicator.isRequestingGradientOf("z_65"):
            communicator.reportGradient("z_65", CoordZ.Gradient(current_design, 65))

# =======================================================================================================
# Perform optimization
# =======================================================================================================

with open("parameters.json",'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

model = KM.Model()

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

    n_iter = int(last_line[0].strip())
    f = float(last_line[1].strip())
    df = float(last_line[2].strip())
    c1 = float(last_line[4].strip())
    c2 = float(last_line[6].strip())
    c3 = float(last_line[8].strip())

    # Check against specifications
    TestCase().assertEqual(n_iter, 10)
    TestCase().assertAlmostEqual(f, 6.61655E+01, 4)
    TestCase().assertAlmostEqual(c1, -9.98583E-01, 4)
    TestCase().assertAlmostEqual(c2, 2.03941E+00, 4)
    TestCase().assertAlmostEqual(c3, 1.00106E+00, 4)

# =======================================================================================================
