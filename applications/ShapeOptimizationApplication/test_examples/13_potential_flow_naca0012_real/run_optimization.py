# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import Kratos core and apps
import KratosMultiphysics
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_response import CreateResponseFunction
import KratosMultiphysics.kratos_utilities as kratos_utils

# Additional imports
from analyzer_base import AnalyzerBaseClass

# Python Libraries
from shutil import copyfile, rmtree
from decimal import *
import math

# Read parameters
with open("optimization_parameters.json",'r') as parameter_file:
    parameters_optimization = Parameters(parameter_file.read())

model = Model()

def solve_potential(current_design, optimization_iteration):
    """ A new directory for the execution of the potential flow analysis is created.
        The *.mdpa file is updated using the updated coordinates from the optimization model part.
        The potential flow analysis is solved in a clean run.
    """

    newpath = "sensitivity_analysis_" + str(optimization_iteration)
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    copyfile("ProjectParametersPrimal.json",  newpath + "/ProjectParametersPrimal.json")
    copyfile("ProjectParametersAdjoint.json",  newpath + "/ProjectParametersAdjoint.json")
    copyfile("response_parameters.json",  newpath + "/response_parameters.json")
    mdpa_name = "naca0012_Case_0_DS_100_AOA_5.0_Far_Field_Mesh_Size_2_Airfoil_Mesh_Size_0.001"
    copyfile(mdpa_name+".mdpa",  newpath + "/"+mdpa_name+".mdpa")

    os.chdir(newpath)
    ## updating the coordinates of the mdpa file

    new_file_name = mdpa_name+".mdpa"
    old_file_name = mdpa_name+"_old.mdpa"
    os.rename(new_file_name, old_file_name)

    if not os.path.isfile(old_file_name):
        raise RuntimeError("mdpa file does not exist", old_file_name)

    new_file = open(new_file_name,'w')
    old_file = open(old_file_name,'r')

    tmp_design_nodes = {} #loop with pybind is slow

    full_design_model = model.GetModelPart(parameters_optimization["optimization_settings"]["model_settings"]["model_part_name"].GetString())

    for node in full_design_model.Nodes:
        tmp_design_nodes[node.Id] = [node.X0, node.Y0, node.Z0]

    inNodes = False

    i = 0
    for line in old_file:
        i += 1

        if line.startswith("End Nodes"):
            inNodes = False

        if inNodes is True:
            line_split = line.split()
            node_id = int(line_split[0])

            # coordinate update
            node_x = format(tmp_design_nodes[node_id][0], '.9f')
            node_y = format(tmp_design_nodes[node_id][1], '.9f')
            node_z = format(tmp_design_nodes[node_id][2], '.9f')
            line_ls = [str(node_id), str(node_x), str(node_y), str(node_z)]
            line_to_write = str(3*" ").join(line_ls)
            new_file.write(line_to_write + '\n')
        else:
            new_file.write(line)

        if line.startswith("Begin Nodes"):
            inNodes = True

    new_file.close()
    old_file.close()

    #########
    print("#"*80)
    print("#### Start sensitivity analysis")
    ## run sensitivity analysis
    with open("response_parameters.json",'r') as parameter_file:
        parameters_sensitivity = KratosMultiphysics.Parameters( parameter_file.read())

    model_sensitivity = KratosMultiphysics.Model()
    response_function = CreateResponseFunction("dummy", parameters_sensitivity["kratos_response_settings"], model_sensitivity)

    response_function.RunCalculation(calculate_gradient=True)
    print("#### End sensitivity analysis")
    print("#"*80)

    ## deleting hdf5 and .res files
    model_part_name = "primal_output"
    for name in os.listdir():
        if name.find(model_part_name) == 0:
            kratos_utils.DeleteFileIfExisting(name)

    os.chdir("..")

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
