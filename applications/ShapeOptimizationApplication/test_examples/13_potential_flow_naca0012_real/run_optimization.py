# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import Kratos core and apps
import KratosMultiphysics as km
import KratosMultiphysics.ShapeOptimizationApplication as kso
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_response import CreateResponseFunction
from KratosMultiphysics.KratosUnittest import WorkFolderScope

# Additional imports
from analyzer_base import AnalyzerBaseClass

# Python Libraries
from shutil import copyfile, rmtree
import os
import decimal
import math

# Read parameters
with open("optimization_parameters.json",'r') as parameter_file:
    parameters_optimization = km.Parameters(parameter_file.read())

model = km.Model()

def update_mdpa_coordinates(mdpa_file_path, current_model_part):
    new_file_name = mdpa_file_path
    old_file_name = mdpa_file_path + ".tmp"
    os.rename(new_file_name, old_file_name)

    with open(new_file_name,'w') as new_file, open(old_file_name,'r') as old_file:

        in_nodes = False

        for i, line in enumerate(old_file):

            if line.startswith("End Nodes"):
                in_nodes = False

            if in_nodes is True:
                line_split = line.split()
                node_id = int(line_split[0])

                # coordinate update
                updated_node = current_model_part.GetNode(node_id)
                node_x = format(updated_node.X0, '.9f')
                node_y = format(updated_node.Y0, '.9f')
                node_z = format(updated_node.Z0, '.9f')
                line_ls = [str(node_id), str(node_x), str(node_y), str(node_z)]
                line_to_write = str(3*" ").join(line_ls)
                new_file.write(line_to_write + '\n')
            else:
                new_file.write(line)

            if line.startswith("Begin Nodes"):
                in_nodes = True

def solve_potential(current_design, optimization_iteration):
    """ A new directory for the execution of the potential flow analysis is created.
        The *.mdpa file is updated using the updated coordinates from the optimization model part.
        The potential flow analysis is solved in a clean run.

        This has some overhead because of the reimport of the model part all the time,
        but it is closer to a standard analysis.
    """

    newpath = "sensitivity_analysis_" + str(optimization_iteration)
    if os.path.exists(newpath):
        shutil.rmtree(newpath)
    os.mkdir(newpath)

    copyfile("ProjectParametersPrimal.json",  os.path.join(newpath, "ProjectParametersPrimal.json"))
    copyfile("ProjectParametersAdjoint.json",  os.path.join(newpath, "ProjectParametersAdjoint.json"))
    copyfile("response_parameters.json",  os.path.join(newpath, "response_parameters.json"))
    mdpa_file_name = "naca0012_Case_0_DS_100_AOA_5.0_Far_Field_Mesh_Size_2_Airfoil_Mesh_Size_0.001" + ".mdpa"
    copyfile(mdpa_file_name,  os.path.join(newpath, mdpa_file_name))

    with WorkFolderScope(newpath, __file__):

        current_model_part = model.GetModelPart(parameters_optimization["optimization_settings"]["model_settings"]["model_part_name"].GetString())
        update_mdpa_coordinates(mdpa_file_name, current_model_part)

        print("#"*80)
        print("Start sensitivity analysis\n")

        ## run sensitivity analysis
        with open("response_parameters.json",'r') as parameter_file:
            parameters_sensitivity = km.Parameters( parameter_file.read())

        model_sensitivity = km.Model()
        response_function = CreateResponseFunction("dummy", parameters_sensitivity["kratos_response_settings"], model_sensitivity)

        response_function.RunCalculation(calculate_gradient=True)

        print("\nEnd sensitivity analysis")
        print("#"*80)

    return response_function.GetValue(), response_function.GetShapeGradient()

# Definition of external analyzer
class CustomAnalyzer(AnalyzerBaseClass):

    def InitializeBeforeOptimizationLoop(self):
        for name in os.listdir():
            if name.find("sensitivity_analysis_") == 0:
                if os.path.isdir(name):
                    rmtree(name)

    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):

        potential_value, potential_gradient = solve_potential(current_design, optimization_iteration)

        # Get the objective function response value
        if communicator.isRequestingValueOf("potential_lift_name"):
            ## get the value of the strain energy
            communicator.reportValue("potential_lift_name", potential_value)

        # Gets the objective function response gradient
        if communicator.isRequestingGradientOf("potential_lift_name"):
            communicator.reportGradient("potential_lift_name", potential_gradient)

# Create optimizer and perform optimization
import optimizer_factory
optimizer = optimizer_factory.CreateOptimizer(parameters_optimization["optimization_settings"], model, CustomAnalyzer())
optimizer.Optimize()
