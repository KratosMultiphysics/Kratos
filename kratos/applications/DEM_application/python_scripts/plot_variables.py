from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.PfemSolidMechanicsApplication import *
from KratosMultiphysics.MachiningApplication import *
# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

class variable_plotter:
    
    def __init__(self, model_part, list_of_nodes_ids):
        
        self.list_of_nodes = []
        self.files = []
        self.model_part = model_part
        
        for node in model_part.Nodes:
            for id in list_of_nodes_ids:
                if node.Id == id:
                    self.list_of_nodes.append(node)
                    file_writer = open("variables_for_node_" + str(id) + ".txt", 'w');
                    file_writer.write("Time  DISPLACEMENT_X  DISPLACEMENT_Y  DISPLACEMENT_Z  VELOCITY_X  VELOCITY_Y  VELOCITY_Z  ACCELERATION_X  ACCELERATION_Y  ACCELERATION_Z  ANGULAR_VELOCITY_X  ANGULAR_VELOCITY_Y  ANGULAR_VELOCITY_Z   ANGULAR_ACCELERATION_X  ANGULAR_ACCELERATION_Y  ANGULAR_ACCELERATION_Z \n")
                    self.files.append(file_writer)
                    print("The Id "+str(id)+" was found in the model part")
                    break
                    
        if len(self.list_of_nodes) != len(list_of_nodes_ids):
            print("Some nodal ids could not be found in the model part! Stopping")
            stop                         
               
        self.plot_variables(0.0)
         
         

    def plot_variables(self, time): 
    
        i = 0
        for file_writer in self.files:
            node = self.list_of_nodes[i]
            string = str(time) \
            + "  " + str(node.GetSolutionStepValue(DISPLACEMENT_X)) \
            + "  " + str(node.GetSolutionStepValue(DISPLACEMENT_Y)) \
            + "  " + str(node.GetSolutionStepValue(DISPLACEMENT_Z)) \
            + "  " + str(node.GetSolutionStepValue(VELOCITY_X)) \
            + "  " + str(node.GetSolutionStepValue(VELOCITY_Y)) \
            + "  " + str(node.GetSolutionStepValue(VELOCITY_Z)) \
            + "  " + str(node.GetSolutionStepValue(ACCELERATION_X)) \
            + "  " + str(node.GetSolutionStepValue(ACCELERATION_Y)) \
            + "  " + str(node.GetSolutionStepValue(ACCELERATION_Z)) \
            + "  " + str(node.GetSolutionStepValue(ANGULAR_VELOCITY_X)) \
            + "  " + str(node.GetSolutionStepValue(ANGULAR_VELOCITY_Y)) \
            + "  " + str(node.GetSolutionStepValue(ANGULAR_VELOCITY_Z)) \
            + "  " + str(node.GetSolutionStepValue(ANGULAR_ACCELERATION_X)) \
            + "  " + str(node.GetSolutionStepValue(ANGULAR_ACCELERATION_Y)) \
            + "  " + str(node.GetSolutionStepValue(ANGULAR_ACCELERATION_Z)) \
            + "\n"
            
            file_writer.write(string)
            i = i + 1 
            
            
            
    def close_files(self):
        
        for file_writer in self.files:
            file_writer.close()
            
        
                
    
