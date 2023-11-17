import KratosMultiphysics
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
    print('IGA ci siamo')

    name_mdpa_true_boundary = "mdpa_files/Weird_shape3" 
    file_mdpa_exists = os.path.isfile(os.path.join(name_mdpa_true_boundary + ".mdpa"))

    if (file_mdpa_exists) :
        current_model = KratosMultiphysics.Model()
        skin_model_part = current_model.CreateModelPart("skin_model_part")
        skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        KratosMultiphysics.ModelPartIO(name_mdpa_true_boundary).ReadModelPart(skin_model_part)
        
        # Write all the points of the skin boundary in an external file
        directory = "txt_files"
        file_name = os.path.join(directory, "true_points.txt")
        if os.path.exists(file_name):
            os.remove(file_name)
        print('ciaooooo3')
        with open(file_name, 'w') as file:
            for condition in skin_model_part.Conditions :
                file.write(f"{condition.GetNodes()[0].X} {condition.GetNodes()[0].Y}\n")
                file.write(f"{condition.GetNodes()[1].X} {condition.GetNodes()[1].Y}\n")
    else :
        if os.path.exists("txt_file/true_points.txt") :
            os.remove("txt_file/true_points.txt")

    def __init__(self, main_model_part, custom_settings):
        super().__init__(main_model_part, custom_settings)
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    def AddVariables(self):
        super().AddVariables()
        # self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SCALAR_LAGRANGE_MULTIPLIER)

    def Initialize(self):
        super().Initialize()

        
    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        # Set all the control_points -> Fix
        # for node in self.main_model_part.GetNodes() :
        #     node.Fix(KratosMultiphysics.TEMPERATURE)
        #     node.Set(KratosMultiphysics.VISITED, True)
        # print('total number of node in main_model_part = ', len(self.main_model_part.GetNodes()  ))

        free_node_x = []
        free_node_y = []
        fixed_node_x = []
        fixed_node_y = []
        # Set Free the active ones
        # with open('Id_active_control_points.txt', 'r') as file:
        #     lines = file.readlines()
        # for line in lines:
        #     node = self.main_model_part.GetNode(int(line))
        #     node.Free(KratosMultiphysics.TEMPERATURE)
        #     node.Set(KratosMultiphysics.VISITED, False)
        #     free_node_x.append(node.X)
        #     free_node_y.append(node.Y)

        # for node in self.main_model_part.GetNodes() :
        #     fixed_node_x.append(node.X)
        #     fixed_node_y.append(node.Y)
        # plt.scatter(free_node_x, free_node_y, marker='o', color='red', label='Free Nodes')
        # plt.scatter(fixed_node_x, fixed_node_y, marker='x', color='black', label='Fixed Nodes')
        
        count = 0
        # for node in self.main_model_part.GetNodes() :
            # plt.annotate(str(node.Id), (node.X, node.Y), textcoords="offset points", xytext=(0,10), ha='center')
            # if node.Is(KratosMultiphysics.VISITED) :
            #     count = count+1
                # plt.annotate(str(node.Id), (node.X, node.Y), textcoords="offset points", xytext=(0,10), ha='center')
                # plt.scatter(node.X, node.Y, marker='v', color='green')
        
        # print('total number of VISITED (Fixed) node in main_model_part = ', count  )
        # plt.scatter(0.0,0.0)
        # plt.scatter(2.0,2.0)

        # with open('txt_files/true_points.txt', 'r') as file:
        #     lines = file.readlines()
        # x_true_values = []
        # y_true_values = []
        # for line in lines:
        #     parts = line.split()  # Dividi la riga in parti utilizzando lo spazio come separatore
        #     if len(parts) == 2:  # Assicurati che ci siano due valori su ogni riga
        #         x_true_values.append(float(parts[0]))
        #         y_true_values.append(float(parts[1]))
        # for i in range(len(x_true_values)-1):
        #     # Inverto che siamo nel physical
        #     plt.plot([y_true_values[i], y_true_values[i+1]],[x_true_values[i], x_true_values[i+1]], color='green')

        # plt.show()
        
   
    
    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        # for node in self.main_model_part.GetNodes() :
        #     print(node.Id, node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
        # exit()
        

    def Finalize(self):
        super().Finalize()





