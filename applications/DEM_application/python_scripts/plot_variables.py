from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

class variable_plotter:
    
    def __init__(self, model_part, list_of_nodes_ids, benchmark_number):
        
        self.list_of_nodes = []
        self.files = []
        self.model_part = model_part
        for node in model_part.Nodes:
            for id in list_of_nodes_ids:
                if node.Id == id:
                    self.list_of_nodes.append(node)
                    file_writer = open("variables_for_node_" + str(benchmark_number) + ".txt", 'w');
                    file_writer.write("#Time  DISPLACEMENT_X  DISPLACEMENT_Y  DISPLACEMENT_Z  ")
                    file_writer.write("ELASTIC_FORCES_X  ELASTIC_FORCES_Y  ELASTIC_FORCES_Z  ")
                    file_writer.write("TOTAL_FORCES_X  TOTAL_FORCES_Y  TOTAL_FORCES_Z  ")
                    file_writer.write("VELOCITY_X  VELOCITY_Y  VELOCITY_Z  ")
                    file_writer.write("ANGULAR_VELOCITY_X  ANGULAR_VELOCITY_Y  ANGULAR_VELOCITY_Z  ")
                    file_writer.write("PARTICLE_MOMENT_X  PARTICLE_MOMENT_Y  PARTICLE_MOMENT_Z\n")
                    self.files.append(file_writer)
                    #print("The Id " + str(id) + " was found in the model part")
                    break
                    
        if len(self.list_of_nodes) != len(list_of_nodes_ids):
            print("Some nodal ids could not be found in the model part! Stopping")                                  
               
        self.plot_variables(0.0)
         
    def plot_variables(self, time): 
        
        i = 0
        
        for file_writer in self.files:
            node = self.list_of_nodes[i]
            string = str(time) \
            + "  " + str(node.GetSolutionStepValue(DISPLACEMENT_X)) \
            + "  " + str(node.GetSolutionStepValue(DISPLACEMENT_Y)) \
            + "  " + str(node.GetSolutionStepValue(DISPLACEMENT_Z)) \
            + "  " + str(node.GetSolutionStepValue(ELASTIC_FORCES_X)) \
            + "  " + str(node.GetSolutionStepValue(ELASTIC_FORCES_Y)) \
            + "  " + str(node.GetSolutionStepValue(ELASTIC_FORCES_Z)) \
            + "  " + str(node.GetSolutionStepValue(TOTAL_FORCES_X)) \
            + "  " + str(node.GetSolutionStepValue(TOTAL_FORCES_Y)) \
            + "  " + str(node.GetSolutionStepValue(TOTAL_FORCES_Z)) \
            + "  " + str(node.GetSolutionStepValue(VELOCITY_X)) \
            + "  " + str(node.GetSolutionStepValue(VELOCITY_Y)) \
            + "  " + str(node.GetSolutionStepValue(VELOCITY_Z)) \
            + "  " + str(node.GetSolutionStepValue(ANGULAR_VELOCITY_X)) \
            + "  " + str(node.GetSolutionStepValue(ANGULAR_VELOCITY_Y)) \
            + "  " + str(node.GetSolutionStepValue(ANGULAR_VELOCITY_Z)) \
            + "  " + str(node.GetSolutionStepValue(PARTICLE_MOMENT)[0]) \
            + "  " + str(node.GetSolutionStepValue(PARTICLE_MOMENT)[1]) \
            + "  " + str(node.GetSolutionStepValue(PARTICLE_MOMENT)[2]) \
            + '\n'
            file_writer.write(string)
            i = i + 1 
               
    def close_files(self):
        
        for file_writer in self.files:
            file_writer.close()
            
class tangential_force_plotter:
    
    def __init__(self, model_part, list_of_nodes_ids, iteration):
        
        self.list_of_nodes = []
        self.files = []
        self.model_part = model_part
        
        for node in model_part.Nodes:
            for id in list_of_nodes_ids:
                if node.Id == id:
                    self.list_of_nodes.append(node)
                    file_writer = open("variables_for_node_" + str(id) + "_iter_" + str(iteration) + ".txt", 'w');
                    file_writer.write("#Time  TOTAL_FORCES_Y  TOTAL_FORCES_Z  ANGULAR_VELOCITY_X\n")
                    self.files.append(file_writer)
                    print("The Id " + str(id) + " was found in the model part")
                    break
                    
        if len(self.list_of_nodes) != len(list_of_nodes_ids):
            print("Some nodal ids could not be found in the model part! Stopping")
                                     
               
        self.plot_tangential_force(0.0)
         
    def plot_tangential_force(self, time): 
        
        i = 0
        
        for file_writer in self.files:
            node = self.list_of_nodes[i]
            string = str(time) \
            + "  " + str(node.GetSolutionStepValue(TOTAL_FORCES_Y)) \
            + "  " + str(node.GetSolutionStepValue(TOTAL_FORCES_Z)) \
            + "  " + str(node.GetSolutionStepValue(ANGULAR_VELOCITY_X)) \
            + '\n'
            file_writer.write(string)
            i = i + 1 
               
    def close_files(self):
        
        for file_writer in self.files:
            file_writer.close()            
