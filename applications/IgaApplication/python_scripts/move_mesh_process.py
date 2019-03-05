from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import numpy as np
# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.IgaApplication import *

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MoveMeshProcess(Model, settings["Parameters"])

class MoveMeshProcess(KratosMultiphysics.Process):

    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters('''{
            "model_part_name"   : "",
            "3d_curve" : {
              "degree": 0,
              "knot_vector": [ ],
              "active_range": [ ],
              "control_points": [[ 1, [0.0,0.0]]]
              },
            "reference_parameter"       : 1.0,
            "current_parameter"         : 1.0
        }''')

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)
        self.model_part = model[self.params["model_part_name"].GetString()]
        
        degree = self.params["3d_curve"]["degree"].GetInt()
        number_cps = self.params["3d_curve"]["control_points"].size()
        
        
        geometry = CurveGeometry3D(degree, number_cps, True)
        

        knot_vector = self.params["3d_curve"]["knot_vector"].GetVector()

        # for i in range(len(knot_vector)): 
        #     geometry.SetKnot(i, knot_vector[i])
        #     print("Knot {}: {}" .format(i, knot_vector[i]))    
            
        
        # geometry.SetKnot(0, 0.0)
        # geometry.SetKnot(1, 0.0)
        # geometry.SetKnot(2, 0.0)
        # geometry.SetKnot(3, 0.0)
        # geometry.SetKnot(4, 40.0)
        # geometry.SetKnot(5, 40.0)
        # geometry.SetKnot(6, 40.0)
        # geometry.SetKnot(7, 40.0)
        
        self.nodes = np.ones((4,number_cps))

        for i in range(number_cps):
            geometry.SetPole(
            self.params["3d_curve"]["control_points"][i][0].GetInt() - 1,
            [self.params["3d_curve"]["control_points"][i][1].GetVector()[0],
            self.params["3d_curve"]["control_points"][i][1].GetVector()[1],
            self.params["3d_curve"]["control_points"][i][1].GetVector()[2]])
            geometry.SetWeight(
                self.params["3d_curve"]["control_points"][i][0].GetInt() - 1,
                self.params["3d_curve"]["control_points"][i][1].GetVector()[3])
            
            print("Pole {}: {}" .format(i, geometry.Pole(i)))
            
            for j in range(3):
                self.nodes[j][i] = self.params["3d_curve"]["control_points"][i][1].GetVector()[j]
        
        print("nodes: \n{}" .format(self.nodes))

            
            
            

        self.reference_location = geometry.PointAt(
            self.params["reference_parameter"].GetDouble())
        
        self.reference_direction = geometry.DerivativesAt(
            self.params["reference_parameter"].GetDouble(), 1)
        
        self.current_location = geometry.PointAt(
            self.params["current_parameter"].GetDouble())
        
        self.current_direction = geometry.DerivativesAt(
            self.params["current_parameter"].GetDouble(), 1)

        

        print("self.reference_direction: {}" .format(self.reference_direction))
        print("self.reference_location: {}" .format(self.reference_location))
        print("self.current_direction: {}" .format(self.current_direction))
        print("self.current_location: {}" .format(self.current_location))
        
        
        # self.reference_location = [0, 0, 0]
        # self.reference_direction = [0, 0, 0]
        

    #def ExecuteInitialize(self):
    #    self.multiple_points_output_process.ExecuteInitialize()

    #def ExecuteBeforeSolutionLoop(self):
    #    self.multiple_points_output_process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        print("ExecuteInitializeSolutionStep")
        
        # translation = [i - j for i, j in zip (self.reference_location, self.current_location)]
        
        # trans_mat = np.identity(4)
        # for i in range(3):
        #     trans_mat[i,3] = translation[i]
        



        # new_coords = np.matmul(trans_mat, self.nodes)
        #for node in self.model_part.Nodes():
        #    new_coords = node.initial_coords * irgen
        #    node.X() = new_coords[0]


    #def ExecuteFinalizeSolutionStep(self):
    #    self.multiple_points_output_process.ExecuteFinalizeSolutionStep()

    #def ExecuteBeforeOutputStep(self):
    #    self.multiple_points_output_process.ExecuteBeforeOutputStep()

    #def ExecuteAfterOutputStep(self):
    #    self.multiple_points_output_process.ExecuteAfterOutputStep()

    #def ExecuteFinalize(self):
    #    self.multiple_points_output_process.ExecuteFinalize()

