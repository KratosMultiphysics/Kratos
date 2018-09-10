from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ContactMechanicsApplication as KratosContact

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

import sys

def CreateParametricWall(main_model_part, custom_settings):
    return ParametricWall(main_model_part, custom_settings)

class ParametricWall(object):

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
            "python_module": "parametric_wall",
            "model_part_name" : "WallDomain",
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
            },
            "bounding_box_settings":{
               "kratos_module": "KratosMultiphysics.ContactMechanicsApplication",
               "bounding_box_type": "SpatialBoundingBox",
               "bounding_box_parameters":{
                   "parameters_list":[],
                   "velocity" : [0.0, 0.0, 0.0],
                   "plane_size": 1.0
               }
            },
            "contact_search_settings":{
               "kratos_module": "KratosMultiphysics.ContactMechanicsApplication",
               "contatct_search_type": "ParametricWallContactSearch",
               "contact_parameters":{
                   "contact_condition_type": "PointContactCondition2D1N",
                   "friction_law_type": "FrictionLaw",
                   "variables_of_properties":{
                     "FRICTION_ACTIVE": false,
                     "MU_STATIC": 0.3,
                     "MU_DYNAMIC": 0.2,
                     "PENALTY_PARAMETER": 1000,
                     "TANGENTIAL_PENALTY_RATIO": 0.1,
                     "TAU_STAB": 1
                   }
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

    def BuildParametricWall(self):

        # construct rigid wall // it will contain the array of nodes, array of elements, and the array of conditions
        self.wall_model_part = self.main_model_part.GetSubModelPart(self.settings["model_part_name"].GetString())
        self.wall_model_part.Set(KratosMultiphysics.RIGID)
        for node in self.wall_model_part.Nodes:
            node.Set(KratosMultiphysics.RIGID,True)
            node.Set(KratosMultiphysics.BOUNDARY,False)

        for node in self.wall_model_part.Conditions:
            node.Set(KratosMultiphysics.ACTIVE,False)

        box_module    = self.settings["bounding_box_settings"]["kratos_module"].GetString()
        box_type_name = self.settings["bounding_box_settings"]["bounding_box_type"].GetString()

        #import module if not previously imported
        module = __import__(box_module)
        module_name = (box_module.split("."))[-1]
        BoundingBox = getattr(getattr(module, module_name), box_type_name)

        #check for the bounding box of a compound wall
        box_settings = self.settings["bounding_box_settings"]["bounding_box_parameters"]

        self.wall_bounding_box = BoundingBox(self.settings["bounding_box_settings"]["bounding_box_parameters"])

        if( self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False ):

            # contruct parametric wall mesh
            self.CreateBoundingBoxMesh(self.wall_bounding_box, self.wall_model_part)

            # construct rigid element // must pass an array of nodes to the element, create a node (CG) and a rigid element set them in the model_part, set the node CG as the reference node of the wall_bounding_box, BLOCKED, set in the wall_model_part for imposed movements processes.
            creation_utility = KratosContact.RigidBodyCreationUtility()
            creation_utility.CreateRigidBodyElement(self.main_model_part, self.wall_bounding_box, self.settings["rigid_body_settings"])

            # create a contact model part
            self.contact_model_part_name =  "contact_"+self.settings["model_part_name"].GetString()

            #can not be a child of wall_model_part due to process imposed variables
            self.main_model_part.CreateSubModelPart(self.contact_model_part_name)
            self.contact_wall_model_part = self.main_model_part.GetSubModelPart(self.contact_model_part_name)

        else:

            # next must be tested:
            self.rigid_body_model_part_name = self.settings["rigid_body_settings"]["rigid_body_model_part_name"].GetString()
            self.rigid_body_model_part = self.main_model_part.GetSubModelPart(self.rigid_body_model_part_name)

            #RigidBodyCenter = self.rigid_body_model_part.GetNode(self.rigid_body_model_part.NumberOfNodes()-1)
            for node in self.rigid_body_model_part.GetNodes():
                RigidBodyCenter = node

            self.wall_bounding_box.SetRigidBodyCenter(RigidBodyCenter);

            # get contact model part
            self.contact_model_part_name =  "contact_"+self.settings["model_part_name"].GetString()

            #can not be a child of wall_model_part due to process imposed variables
            self.contact_wall_model_part = self.main_model_part.GetSubModelPart(self.contact_model_part_name)



        #construct the search strategy
        search_module    = self.settings["contact_search_settings"]["kratos_module"].GetString()
        search_type_name = self.settings["contact_search_settings"]["contact_search_type"].GetString()

        #import module if not previously imported
        smodule = __import__(search_module)
        smodule_name = (search_module.split("."))[-1]
        SearchProcess = getattr(getattr(smodule, smodule_name), search_type_name)

        print("::[Parametric_Wall]:: Contact Model Part",self.contact_model_part_name)

        self.SearchStrategy = SearchProcess(self.main_model_part, self.contact_model_part_name, self.wall_bounding_box, self.settings["contact_search_settings"]["contact_parameters"])

        print("::[Parametric_Wall]:: -BUILT-")

    ####

    #
    def CreateBoundingBoxMesh(self, bounding_box, model_part):

        # construct bounding box mesh of surface elements
        # Note: new nodes must be inserted in boundary conditions subdomains
        number_of_linear_partitions  = 10
        number_of_angular_partitions = 20

        bounding_box.CreateBoundingBoxBoundaryMesh(model_part, number_of_linear_partitions, number_of_angular_partitions)

        # set mesh upper and lower points
        upper_point = KratosMultiphysics.Array3()
        upper = self.GetUpperPoint(model_part)
        print("upper",upper)
        for i in range(0,len(upper)):
            upper_point[i] = upper[i]
            bounding_box.SetUpperPoint(upper_point)

        lower_point = KratosMultiphysics.Array3()
        lower = self.GetLowerPoint(model_part)
        print("lower",lower)
        for i in range(0,len(lower)):
            lower_point[i] = lower[i]
            bounding_box.SetLowerPoint(lower_point)

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
            if( node.Z < min_z ):
                min_z = node.Z

        if( dimension == 2 ):
            return [min_x, min_y, 0]
        else:
            return [min_x, min_y, min_z]

    ####

    def Initialize(self):

        self.SearchStrategy.ExecuteInitialize()


    def InitializeSearch(self):
        pass


    def FinalizeSearch(self):
        pass


    def ExecuteSearch(self):

        self.SearchStrategy.Execute()
