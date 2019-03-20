from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import numpy as np
# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.IgaApplication import *

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MoveObjectProcess(Model, settings["Parameters"])

class MoveObjectProcess(KratosMultiphysics.Process):

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
            "reference_parameter"       : 0.0,
            "steps"                     : 0 
        }''')

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)
        self.model_part = model[self.params["model_part_name"].GetString()]
        
        degree = self.params["3d_curve"]["degree"].GetInt()
        number_cps = self.params["3d_curve"]["control_points"].size()
        
        self.geometry = CurveGeometry3D(degree, number_cps, True)
        
        for i in range(number_cps):
            
            cp_idx = self.params["3d_curve"]["control_points"][i][0].GetInt() - 1
            cp = self.params["3d_curve"]["control_points"][i][1].GetVector()

            self.geometry.SetPole(cp_idx, [cp[0], cp[1], cp[2]])
            self.geometry.SetWeight(cp_idx, cp[3])

        knot_vector = self.params["3d_curve"]["knot_vector"].GetVector()

        for i in range(1,len(knot_vector) - 1): 
            self.geometry.SetKnot(i-1, knot_vector[i])
    
    def ExecuteInitializeSolutionStep(self):
        print("ExecuteInitializeSolutionStep")
        
        if self.params["steps"].GetInt() != 0:
            step_size = ((self.params["reference_parameter"].GetDouble() - self.params["current_parameter"].GetDouble())
                        /(self.params["steps"].GetInt()))
        
        
        index = 0
        coords = np.ones((4, len(self.model_part.Nodes)))   
        for node in self.model_part.Nodes:
            coords[0,index] = node.X0
            coords[1,index] = node.Y0
            coords[2,index] = node.Z0
            index += 1

        current_val = self.geometry.DerivativesAt(self.params["current_parameter"].GetDouble(), 1)

        import matplotlib.pyplot as plt
        plt.plot([coords[0,0],coords[0,1],coords[0,3],coords[0,2],coords[0,0]],
                [coords[1,0],coords[1,1],coords[1,3],coords[1,2],coords[1,0]])
        
        for j in range(4):
            plt.annotate(str(j), xy=(coords[0,j], coords[1,j]),color="red")

        for i in range(self.params["steps"].GetInt()):
            
            next_parameter = (i + 1) * step_size
            next_val = self.geometry.DerivativesAt(next_parameter, 1)
            
            translation_current = current_val[0]
            translation_reference = next_val[0]

            angle = (np.arctan2(next_val[1][1], next_val[1][0])
                    -np.arctan2(current_val[1][1], current_val[1][0]))

            M = np.identity(4)
            M[0,0] =  np.cos(angle)
            M[0,1] = -np.sin(angle)
            M[1,0] =  np.sin(angle)
            M[1,1] =  np.cos(angle)

            M[0,3] = (translation_reference[0] - np.cos(angle) * translation_current[0]
                    + np.sin(angle) * translation_current[1])

            M[1,3] = (translation_reference[1] - np.cos(angle) * translation_current[1] 
                    - np.sin(angle) * translation_current[0])
            
            M[2,3] = translation_reference[2] - translation_current[2]

            coords = np.matmul(M, coords)

            plt.plot([coords[0,0],coords[0,1],coords[0,3],coords[0,2],coords[0,0]],
                    [coords[1,0],coords[1,1],coords[1,3],coords[1,2],coords[1,0]])
            
            for j in range(4):
                plt.annotate(str(j), xy=(coords[0,j], coords[1,j]),color="red")


            current_val = next_val
       
        plt.grid()
        plt.axis("equal")
        plt.show()
