# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import Kratos core and apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *
import KratosMultiphysics.StructuralMechanicsApplication
import structural_response_function_factory

# Additional imports
from analyzer_base import AnalyzerBaseClass

# Python Libraries
from shutil import copyfile
from decimal import *
import math


# Read parameters
with open("optimization_parameters.json",'r') as parameter_file:
    parameters_optimization = Parameters(parameter_file.read())

# Definition of external analyzer
class CustomAnalyzer(AnalyzerBaseClass):

    def InitializeBeforeOptimizationLoop(self):
        pass
        #TODO delete all sensitivity_analysis folders module glob

    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):

        newpath = "sensitivity_analysis_" + str(optimization_iteration)
        if not os.path.exists(newpath):
            os.makedirs(newpath)

        copyfile("plate_analysis_parameters.json",  newpath + "/plate_analysis_parameters.json")
        #copyfile("plate_adjoint_analysis_parameters.json",  newpath + "/plate_adjoint_analysis_parameters.json")
        copyfile("plate_adjoint_analysis_parameters.json",  newpath + "/plate_adjoint_analysis_parameters.json")
        copyfile("adjoint_strain_energy_response_parameters_plate.json",  newpath + "/adjoint_strain_energy_response_parameters_plate.json")
        copyfile("2D_material.json",  newpath + "/2D_material.json")
        copyfile("plate_surface_load.mdpa",  newpath + "/plate_surface_load.mdpa")

        os.chdir(newpath)
        ## updating the coordinates of the mdpa file

        new_file_name = "plate_surface_load.mdpa"
        old_file_name = "plate_surface_load_old.mdpa"
        os.rename(new_file_name, old_file_name)

        if not os.path.isfile(old_file_name):
            raise RuntimeError("mdpa file does not exist", old_file_name)

        new_file = open(new_file_name,'w')
        old_file = open(old_file_name,'r')

        tmp_design_nodes = {} #loop with pybind is slow

        for node in current_design.Nodes:
            tmp_design_nodes[node.Id] = [node.X0, node.Y0, node.Z0]
            if node.Id == 2256:
                print(node)
            if math.isnan(node.X0):
                raise RuntimeError(node.Id, node.X0)
            if math.isnan(node.Y0):
                raise RuntimeError(node.Id, node.Y0)
            if math.isnan(node.Z0):
                raise RuntimeError(node.Id, node.Z0)

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
        with open("adjoint_strain_energy_response_parameters_plate.json",'r') as parameter_file:
            parameters_sensitivity = KratosMultiphysics.Parameters( parameter_file.read())

        model_sensitivity = KratosMultiphysics.Model()
        response_function = structural_response_function_factory.CreateResponseFunction("dummy", parameters_sensitivity["kratos_response_settings"], model_sensitivity)

        model_part_primal = response_function.primal_model_part
        model_part_adjoint = response_function.adjoint_model_part

        response_function.RunCalculation(calculate_gradient=True)
        print("#### End sensitivity analysis")
        print("#"*80)
        os.chdir("..")

        # Get the objective function response value
        if communicator.isRequestingValueOf("strain_energy"):
            ## get the value of the strain energy
            value = response_function.GetValue()
            communicator.reportValue("strain_energy", value)

        # Gets the objective function response gradient
        if communicator.isRequestingGradientOf("strain_energy"):
            gradient = response_function.GetShapeGradient()
            communicator.reportGradient("strain_energy", gradient)

model = Model()
# Create optimizer and perform optimization
import optimizer_factory
optimizer = optimizer_factory.CreateOptimizer(parameters_optimization["optimization_settings"], model, CustomAnalyzer())
optimizer.Optimize()
