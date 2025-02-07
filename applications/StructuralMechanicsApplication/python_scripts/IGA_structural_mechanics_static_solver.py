# Importing the Kratos Library
import KratosMultiphysics

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver
import json
import os
import matplotlib.pyplot as plt
import numpy as np

def CreateSolver(model, custom_settings):
    return IgaStructuralMechanicsStaticSolver(model, custom_settings)

class IgaStructuralMechanicsStaticSolver(MechanicalSolver):
    """The structural mechanics static solver.

    This class creates the mechanical solvers for static analysis.

    See structural_mechanics_solver.py for more information.
    """
    with open('txt_files/input_data.txt', 'r') as file_in:
        line = file_in.readline().strip()
        name_mdpa_true_boundary = line
        # name_mdpa_true_boundary = "mdpa_files/Weird_shape3" 
        file_mdpa_exists = os.path.isfile(os.path.join(name_mdpa_true_boundary + ".mdpa"))
        
        if (file_mdpa_exists) :
            print('mdpa file of embedded boundary exists!')
            current_model = KratosMultiphysics.Model()
            skin_model_part2 = current_model.CreateModelPart("skin_model_part2")
            skin_model_part2.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
            KratosMultiphysics.ModelPartIO(name_mdpa_true_boundary).ReadModelPart(skin_model_part2)
            
            # Write all the points of the skin boundary in an external file
            directory = "txt_files"
            file_name = os.path.join(directory, "true_points.txt")
            if os.path.exists(file_name):
                os.remove(file_name)
            with open(file_name, 'w') as file:
                for condition in skin_model_part2.Conditions :
                    file.write(f"{condition.GetNodes()[0].X} {condition.GetNodes()[0].Y}\n")
                    file.write(f"{condition.GetNodes()[1].X} {condition.GetNodes()[1].Y}\n")
                    
                    
                    
        # REMOVE
        line = file_in.readline().strip()
        line = file_in.readline().strip()
        line = file_in.readline().strip()
        name_mdpa_true_boundary = line
        # name_mdpa_true_boundary = "mdpa_files/Weird_shape3" 
        file_mdpa_exists = os.path.isfile(os.path.join(name_mdpa_true_boundary + ".mdpa"))
        
        if (file_mdpa_exists) :
            print('mdpa file of embedded boundary exists!')
            current_model = KratosMultiphysics.Model()
            skin_model_part2 = current_model.CreateModelPart("skin_model_part2")
            skin_model_part2.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
            KratosMultiphysics.ModelPartIO(name_mdpa_true_boundary).ReadModelPart(skin_model_part2)
            
            # Write all the points of the skin boundary in an external file
            directory = "txt_files"
            file_name = os.path.join(directory, "true_points2.txt")
            if os.path.exists(file_name):
                os.remove(file_name)
            with open(file_name, 'w') as file:
                for condition in skin_model_part2.Conditions :
                    file.write(f"{condition.GetNodes()[0].X} {condition.GetNodes()[0].Y}\n")
                    file.write(f"{condition.GetNodes()[1].X} {condition.GetNodes()[1].Y}\n")
        
        line = file_in.readline().strip()
        name_mdpa_true_boundary = line
        # name_mdpa_true_boundary = "mdpa_files/Weird_shape3" 
        file_mdpa_exists = os.path.isfile(os.path.join(name_mdpa_true_boundary + ".mdpa"))
        
        if (file_mdpa_exists) :
            print('mdpa file of embedded boundary exists!')
            current_model = KratosMultiphysics.Model()
            skin_model_part2 = current_model.CreateModelPart("skin_model_part2")
            skin_model_part2.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
            KratosMultiphysics.ModelPartIO(name_mdpa_true_boundary).ReadModelPart(skin_model_part2)
            
            # Write all the points of the skin boundary in an external file
            directory = "txt_files"
            file_name = os.path.join(directory, "true_points3.txt")
            if os.path.exists(file_name):
                os.remove(file_name)
            with open(file_name, 'w') as file:
                for condition in skin_model_part2.Conditions :
                    file.write(f"{condition.GetNodes()[0].X} {condition.GetNodes()[0].Y}\n")
                    file.write(f"{condition.GetNodes()[1].X} {condition.GetNodes()[1].Y}\n")
                
    
    
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[StaticMechanicalSolver]:: ", "Construction finished")
    
    
    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "Initializing ...")
        
        # The mechanical solution strategy is created here if it does not already exist.
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        mechanical_solution_strategy = self._GetSolutionStrategy()
        print(mechanical_solution_strategy)
        mechanical_solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        mechanical_solution_strategy.Initialize()

        # Printing that inialization is finished
        KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "Finished initialization.")
    
    def FinalizeSolutionStep(self):
        # self.printDofsAndCPs()
        return super().FinalizeSolutionStep()
        
        
    def _GetScheme(self):
        if not hasattr(self, '_solution_scheme'):
            self._solution_scheme = self._CreateScheme()
        return self._solution_scheme

    def _CreateScheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
    



    def printDofsAndCPs(self) :
        fig, ax = plt.subplots()

        free_node_x = []
        free_node_y = []
        free_node_z = []
        fixed_node_x = []
        fixed_node_y = []
        fixed_node_z = []
        dof = []

        # z_ref = 1.0
        # Set Free the active ones
        if os.path.exists("txt_files/Id_active_control_points.txt"):
            with open('txt_files/Id_active_control_points.txt', 'r') as file:
                lines = file.readlines()
            for line in lines:
                numbers = line.split()
                node = self.main_model_part.GetNode(int(numbers[0]))

                free_node_x.append(node.X)
                free_node_y.append(node.Y)
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

                free_node_x2.append(node.X)
                free_node_y2.append(node.Y)
                dof2.append(numbers[1])
        
        for node in self.main_model_part.GetNodes() :
            fixed_node_x.append(node.X)
            fixed_node_y.append(node.Y)
            fixed_node_z.append(node.Z)

        
        # ax = fig.add_subplot(111, projection='3d')
        print(free_node_x)
        print(free_node_y)

        # ax.scatter(fixed_node_x, fixed_node_y, fixed_node_z, marker='x', color='black', s=12,  label='Fixed Nodes')
        ax.scatter(free_node_x, free_node_y, marker='d', color='red',s=60, label='Free Nodes')

        ax.scatter(free_node_x2, free_node_y2, marker='x', alpha=0.3, color='green',s=25, label='Free Nodes')
        
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

        # plt.show()
