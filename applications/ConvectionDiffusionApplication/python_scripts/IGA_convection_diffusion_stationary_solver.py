import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import numpy as np
# from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import json
import os

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_stationary_solver

def CreateSolver(main_model_part, custom_settings):
    return IGAConvectionDiffusionStationarySolver(main_model_part, custom_settings)


class IGAConvectionDiffusionStationarySolver(convection_diffusion_stationary_solver.ConvectionDiffusionStationarySolver):

    # Write all the points of the skin boundary in an external file
    directory = "txt_files"
    file_name = os.path.join(directory, "true_points.txt")
    if os.path.exists(file_name):
        os.remove(file_name)
    name_mdpa_true_boundary = "mdpa_files/Weird_shape3"
    # name_mdpa_true_boundary = "mdpa_files/sphere" 
    file_mdpa_exists = os.path.isfile(os.path.join(name_mdpa_true_boundary + ".mdpa"))
    
    if (file_mdpa_exists) :
        print('mdpa file of embedded boundary exists!')
        current_model = KratosMultiphysics.Model()
        skin_model_part2 = current_model.CreateModelPart("skin_model_part2")
        skin_model_part2.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        KratosMultiphysics.ModelPartIO(name_mdpa_true_boundary).ReadModelPart(skin_model_part2)

        with open(file_name, 'w') as file:
            for condition in skin_model_part2.Conditions :
                if (condition.GetGeometry().WorkingSpaceDimension() <= 2) :
                    file.write(f"{condition.GetNodes()[0].X} {condition.GetNodes()[0].Y}\n")
                    file.write(f"{condition.GetNodes()[1].X} {condition.GetNodes()[1].Y}\n")
                else :
                    file.write(f"{condition.GetNodes()[0].X} {condition.GetNodes()[0].Y} {condition.GetNodes()[0].Z} \n")
                    file.write(f"{condition.GetNodes()[1].X} {condition.GetNodes()[1].Y} {condition.GetNodes()[1].Z} \n")
                    file.write(f"{condition.GetNodes()[2].X} {condition.GetNodes()[2].Y} {condition.GetNodes()[2].Z} \n")
                    
    

    # Write all the points of the skin boundary in an external file
    directory = "txt_files"
    file_name = os.path.join(directory, "true_points_outer.txt")
    if os.path.exists(file_name):
        os.remove(file_name)
    name_mdpa_true_boundary = "mdpa_files/external_bunny.."
    file_mdpa_exists = os.path.isfile(os.path.join(name_mdpa_true_boundary + ".mdpa"))
    if (file_mdpa_exists) :
        current_model = KratosMultiphysics.Model()
        skin_model_part2_outer = current_model.CreateModelPart("skin_model_part2_outer")
        skin_model_part2_outer.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        KratosMultiphysics.ModelPartIO(name_mdpa_true_boundary).ReadModelPart(skin_model_part2_outer)
        
        with open(file_name, 'w') as file:
            for condition in skin_model_part2_outer.Conditions :
                file.write(f"{condition.GetNodes()[0].X} {condition.GetNodes()[0].Y}\n")
                file.write(f"{condition.GetNodes()[1].X} {condition.GetNodes()[1].Y}\n")


    
    def __init__(self, main_model_part, custom_settings):
        super().__init__(main_model_part, custom_settings)
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    def AddVariables(self):
        super().AddVariables()
        # self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_X)
        # self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_Z)
        # self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_Y)
        # self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SCALAR_LAGRANGE_MULTIPLIER)

    def AddDofs(self):
        super().AddDofs()
        # KratosMultiphysics.VariableUtils().AddDof(KratosConvDiff.SCALAR_LAGRANGE_MULTIPLIER, self.main_model_part)


    def Initialize(self):
        super().Initialize()

        
    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        
    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        # self.printDofsAndCPs() # Print control points and dofs

    def Finalize(self):
        super().Finalize()




















    def printDofsAndCPs(self) :
        # fig, ax = plt.subplots()

        free_node_x = []
        free_node_y = []
        free_node_z = []
        fixed_node_x = []
        fixed_node_y = []
        fixed_node_z = []
        dof = []

        z_ref = 1.0
        # Set Free the active ones
        if os.path.exists("txt_files/Id_active_control_points.txt"):
            with open('txt_files/Id_active_control_points.txt', 'r') as file:
                lines = file.readlines()
            for line in lines:
                numbers = line.split()
                node = self.main_model_part.GetNode(int(numbers[0]))
                node.Free(KratosMultiphysics.TEMPERATURE)
                node.Set(KratosMultiphysics.VISITED, False)


                if (np.abs(node.Y - z_ref) > 1e-1): continue 
                free_node_x.append(node.X)
                free_node_y.append(node.Y)
                free_node_z.append(node.Z)
                dof.append(numbers[1])
        
        dof2 = []
        free_node_x2 = []
        free_node_y2= []
        free_node_z2= []
        # Set Free the active ones
        if os.path.exists("txt_files/Id_active_control_points_condition.txt"):
            with open('txt_files/Id_active_control_points_condition.txt', 'r') as file:
                lines = file.readlines()
            for line in lines:
                numbers = line.split()
                node = self.main_model_part.GetNode(int(numbers[0]))

                if (np.abs(node.Y - z_ref) > 1e-1): continue 
                free_node_x2.append(node.X)
                free_node_y2.append(node.Y)
                free_node_z2.append(node.Z)
                dof2.append(numbers[1])
        
        for node in self.main_model_part.GetNodes() :

            if (np.abs(node.Y - z_ref) > 1e-1): continue 
            fixed_node_x.append(node.X)
            fixed_node_y.append(node.Y)
            fixed_node_z.append(node.Z)

        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # ax.scatter(fixed_node_x, fixed_node_y, fixed_node_z, marker='x', color='black', s=12,  label='Fixed Nodes')
        ax.scatter(free_node_x, free_node_y, free_node_z, marker='d', color='red',s=8, label='Free Nodes')

        ax.scatter(free_node_x2, free_node_y2, free_node_z2, marker='x', alpha=0.3, color='green',s=12, label='Free Nodes')
        
        # count = 0
        # for i in range(len(dof)) :
        #     ax.annotate(str(dof[i]), (free_node_x[i], free_node_y[i], free_node_z[i]), textcoords="offset points", xytext=(0,10), ha='center')

        # for i in range(len(dof2)) :
        #     ax.annotate(str(dof2[i]), (free_node_x2[i], free_node_y2[i], free_node_z2[i]), textcoords="offset points", color='green', xytext=(0,8), ha='center')
     
        # Read the refinements.iga.json
        # with open('refinements.iga.json', 'r') as file:
        #     refinements_parameters = json.load(file)
        # insert_nb_per_span_u = refinements_parameters['refinements'][0]['parameters']['insert_nb_per_span_u']
        # insert_nb_per_span_v = refinements_parameters['refinements'][0]['parameters']['insert_nb_per_span_v']

        # knots_u = [0.0]
        # knots_v = [0.0]
        # initial = 0.0
        # total = 2.0
        # for j in range(1, insert_nb_per_span_u + 1):
        #     knots_u.append(initial + total / (insert_nb_per_span_u + 1) * j)
        # for j in range(1, insert_nb_per_span_v + 1):
        #     knots_v.append(initial + total / (insert_nb_per_span_v + 1) * j)
        # knots_u.append(2.0)
        # knots_v.append(2.0)

        # ax.set_xlim(min(knots_u)-0.2, max(knots_u)+0.2)
        # ax.set_ylim(min(knots_v)-0.2, max(knots_v)+0.2)

        # for u in knots_u:
        #     ax.axvline(x=u, color='red', linestyle='--', linewidth=1.0)
        # for v in knots_v:
        #     ax.axhline(y=v, color='red', linestyle='--', linewidth=1.0)
        # ax.grid(True, linestyle='--', linewidth=0.5)
        # ax.gca().set_aspect('equal', adjustable='box')

        plt.show()


