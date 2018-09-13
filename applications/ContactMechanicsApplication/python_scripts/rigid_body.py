from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.ContactMechanicsApplication as KratosContact

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

import sys

def CreateRigidBody(main_model_part, custom_settings):
    return RigidBody(main_model_part, custom_settings)

class RigidBody(object):

    ##constructor. the constructor shall only take care of storing the settings
    ##and the pointer to the main_model part.
    ##
    ##real construction shall be delayed to the function "Initialize" which
    ##will be called once the mesher is already filled
    def __init__(self, main_model_part, custom_settings):

        self.main_model_part = main_model_part

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "python_module"   : "rigid_body",
            "model_part_name" : "RigidBodyDomain",
            "rigid_body_settings":{
               "rigid_body_element_type": "TranslatoryRigidElement3D1N",
               "fixed_body": true,
               "compute_body_parameters": false,
               "rigid_body_model_part_name": "RigidBodyDomain",
               "rigid_body_parameters":{
                   "center_of_gravity": [0.0 ,0.0, 0.0],
                   "mass":0.0,
                   "main_inertias": [0.0, 0.0, 0.0],
                   "main_axes": [ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ]
               }
            }
        }
        """)

        ## new node and rigid body element inside the same mesh : boundary conditions also applied
        ## this node and elements must be considered in the computing model part
        ## new contact conditions must be already assembled

        ## if exist a movement from a point different from CG a link condition must be used

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        # construct rigid body // it will contain the array of nodes, array of elements, and the array of conditions
        self.rigid_body_model_part = self.main_model_part.GetSubModelPart(self.settings["model_part_name"].GetString())
        self.rigid_body_model_part.Set(KratosMultiphysics.RIGID)

        for node in self.rigid_body_model_part.Nodes:
            node.Set(KratosMultiphysics.RIGID,True)

        for node in self.rigid_body_model_part.Conditions:
            node.Set(KratosMultiphysics.ACTIVE,False)


        #check for the bounding box of a compound wall
        box_settings = KratosMultiphysics.Parameters("""
        {
            "parameters_list":[{
                   "upper_point": [0.0, 0.0, 0.0],
                   "lower_point": [0.0, 0.0, 0.0],
                   "convexity": 1
            }],
            "velocity" : [0.0, 0.0, 0.0]
        }
        """)

        box_parameters = box_settings["parameters_list"][0]

        upper_point = self.GetUpperPoint(self.rigid_body_model_part)
        counter = 0
        for i in upper_point:
            box_parameters["upper_point"][counter].SetDouble(i)
            counter+=1

        lower_point = self.GetLowerPoint(self.rigid_body_model_part)
        counter = 0
        for i in lower_point:
            box_parameters["lower_point"][counter].SetDouble(i)
            counter+=1

        self.bounding_box = KratosDelaunay.SpatialBoundingBox(box_settings)

        # construct rigid element // must pass an array of nodes to the element, create a node (CG) and a rigid element set them in the model_part, set the node CG as the reference node of the wall_bounding_box, BLOCKED, set in the wall_model_part for imposed movements processes.
        creation_utility = KratosContact.RigidBodyCreationUtility()
        creation_utility.CreateRigidBodyElement(self.main_model_part, self.bounding_box, self.settings["rigid_body_settings"])


        print("::[Rigid_Body]:: -BUILT-")

    ####

    #
    def GetUpperPoint(self, model_part):

        dimension = model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION]

        max_x = sys.float_info.min
        max_y = sys.float_info.min
        max_z = sys.float_info.min

        for node in model_part.Nodes:
            if( node.X > max_x ):
                max_x = node.X
            if( node.Y > max_y ):
                max_y = node.Y
            if( node.Z > max_z ):
                max_z = node.Z

        if( dimension == 2 ):
            return [max_x, max_y, 0]
        else:
            return [max_x, max_y, max_z]

    #
    def GetLowerPoint(self, model_part):

        dimension = model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION]

        min_x = sys.float_info.max
        min_y = sys.float_info.max
        min_z = sys.float_info.max

        for node in model_part.Nodes:
            if( node.X < min_x ):
                min_x = node.X
            if( node.Y < min_y ):
                min_y = node.Y
            if( node.Z > min_z ):
                min_z = node.Z

        if( dimension == 2 ):
            return [min_x, min_y, 0]
        else:
            return [min_x, min_y, min_z]

    ####

    def Initialize(self):
        pass
