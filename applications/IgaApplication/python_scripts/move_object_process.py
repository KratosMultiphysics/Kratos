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
            "current_parameter"         : 0.0,
            "reference_parameter"       : 0.0
        }''')

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)
        self.model_part = model[self.params["model_part_name"].GetString()]
        
        degree = self.params["3d_curve"]["degree"].GetInt()
        number_cps = self.params["3d_curve"]["control_points"].size()
        
        geometry = CurveGeometry3D(degree, number_cps, True)
        
        for i in range(number_cps):
            
            cp_idx = self.params["3d_curve"]["control_points"][i][0].GetInt() - 1
            cp = self.params["3d_curve"]["control_points"][i][1].GetVector()

            geometry.SetPole(cp_idx, [cp[0], cp[1], cp[2]])
            geometry.SetWeight(cp_idx, cp[3])

        knot_vector = self.params["3d_curve"]["knot_vector"].GetVector()

        for i in range(1,len(knot_vector) - 1): 
            geometry.SetKnot(i-1, knot_vector[i])
            


        reference_val = geometry.DerivativesAt(
            self.params["reference_parameter"].GetDouble(), 1)
        current_val = geometry.DerivativesAt(
            self.params["current_parameter"].GetDouble(), 1)
        
        self.reference_location = reference_val[0]
        self.reference_direction = reference_val[1]
        self.current_location = current_val[0]
        self.current_direction = current_val[1]
        
    
    def ExecuteInitializeSolutionStep(self):
        print("ExecuteInitializeSolutionStep")
        
        translation_origin = [-i for i in self.current_location]
        translation_reference = self.reference_location
        
        angle = np.arccos(
            np.dot(self.reference_direction, self.current_direction) /
            (np.linalg.norm(self.reference_direction) * np.linalg.norm(self.current_direction)))
        

        coords = np.ones((4, len(self.model_part.Nodes)))        
        i = 0
        for node in self.model_part.Nodes:
            coords[0,i] = node.X0
            coords[1,i] = node.Y0
            coords[2,i] = node.Z0
            i += 1

        M = np.identity(4)
        M[0,0] = np.cos(angle)
        M[0,1] = -np.sin(angle)
        M[1,0] = np.sin(angle)
        M[1,1] = np.cos(angle)
        
        M[0,3] = (translation_origin[0] * np.cos(angle) 
                - np.sin(angle) * translation_origin[1] 
                + translation_reference[0])
                
        M[1,3] = (translation_origin[0] * np.sin(angle) 
                - np.cos(angle) * translation_origin[1] 
                + translation_reference[1])

        M[2,3] = translation_origin[2] + translation_reference[2]

        new_coords = np.matmul(M, coords)


        import matplotlib.pyplot as plt
        

        plt.plot([coords[0,0],coords[0,1],coords[0,3],coords[0,2],coords[0,0]],
                 [coords[1,0],coords[1,1],coords[1,3],coords[1,2],coords[1,0]], 
                 color="black")
        for i in range(4):
            plt.annotate(str(i), xy=(coords[0,i], coords[1,i]),color="black")

        plt.plot([new_coords[0,0],new_coords[0,1],new_coords[0,3],new_coords[0,2],new_coords[0,0]],
                 [new_coords[1,0],new_coords[1,1],new_coords[1,3],new_coords[1,2],new_coords[1,0]], 
                 color="red")
        for i in range(4):
            plt.annotate(str(i), xy=(new_coords[0,i], new_coords[1,i]),color="red")


        print("Current Location: {}" .format(self.current_location))
        print("Current Direction: {}" .format(self.current_direction))
        print("Reference Location: {}" .format(self.reference_location))
        print("Reference Direction: {}" .format(self.reference_direction))

        print("Initial Coordinates; \n{}" .format(coords))
        print("Transformed Coordinates; \n{}" .format(new_coords))
        plt.axis([-40, 40, -40, 40])
        plt.show()


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

