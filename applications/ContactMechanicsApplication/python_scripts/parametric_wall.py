from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ContactMechanicsApplication as KratosContact

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateParametricWall(Model, main_model_part, custom_settings):
    return ParametricWall(Model, main_model_part, custom_settings)

class ParametricWall(object):
    
    ##constructor. the constructor shall only take care of storing the settings 
    ##and the pointer to the main_model part.
    ##
    ##real construction shall be delayed to the function "Initialize" which 
    ##will be called once the modeler is already filled
    def __init__(self, Model, model_part, custom_settings):
        
        self.model_part = model_part    
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "python_file_name": "parametric_wall",
            "mesh_id": 0,
            "sub_model_part_name" : "Wall Domain",
            "rigid_body_settings":{
               "rigid_body_element_type": "TranslatoryRigidElement3D1N",
               "fixed_body": true,
               "compute_parameters": false,
               "rigid_body_parameters":{
                   "center_of_gravity": [0.0 ,0.0, 0.0],
                   "mass":0.0,
                   "main_inertias": [0.0, 0.0, 0.0],
                   "main_axes": [ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ]
               }
            },
            "bounding_box_settings":{
               "implemented_in_module": "KratosMultiphysics.ContactMechanisApplication",
               "bounding_box_type": "SpatialBoundingBox",
               "bounding_box_parameters":{
                   "parameters_list":[],
                   "velocity" : [0.0, 0.0, 0.0]
               }
            },
            "contact_search_settings":{
               "python_file_name": "parametric_wall_contact_search",
               "search_frequency": 0,            
               "contact_parameters":{
                   "contact_condition_type": "PointContactCondition2D1N",
                   "friction_law_type": "FrictionLaw",
                   "implemented_in_module": "KratosMultiphysics.ContactMechanicsApplication",
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
                
        #construct rigid wall // it will contain the array of nodes, array of elements, and the array of conditions
        self.wall_model_part = Model[self.settings["sub_model_part_name"].GetString()]

        module = __import__(self.settings["bounding_box_settings"]["implemented_in_module"].GetString())      
        box_module    = self.settings["bounding_box_settings"]["implemented_in_module"].GetString()
        box_type_name = self.settings["bounding_box_settings"]["bounding_box_type"].GetString()
        BoundingBox   = None

        box_type_call = " BoundingBox = " box_module + box_type_name
        exec(box_type_call)
        
        self.wall_bounding_box = BoundingBox(self.wall_bounding_box, self.settings["bounding_box_settings"]["bounding_box_parameters"])

        #construct rigid element // must pass an array of nodes to the element, create a node (CG) and a rigid element set them in the model_part, set the node CG as the reference node of the wall_bounding_box, BLOCKED, set in the wall_model_part for imposed movements processes.
        creation_utility = KratosContact.RigidBodyCreationUtilities()
        creation_utility.CreateRigidBodyElement(self.model_part, self.wall_bounding_box, self.settings["rigid_body_settings"])

        #construct the search strategy
        search_module = __import__(self.settings["contact_search_strategy"]["python_file_name"].GetString())
        self.SearchStrategy = search_module.CreateContactSearch(self.model_part, self.wall_bounding_box, self.settings["contact_search_settings"])

        print("Construction of Parametric Wall finished")
        

    #### 

    def Initialize(self):

        print("::[Parametric_Wall]:: -START-")
        
        self.domain_size = self.model_part.ProcessInfo[DOMAIN_SIZE]
        self.mesh_id     = self.settings["mesh_id"].GetInt()

        # Meshing Stratety
        self.SearchStrategy.Initialize()   #in the contact search a contact_mesh or model_part must be created for each parametric wall.

        # create next inside the SearchStrategy: c++
        self.contact_model_part_name =  "contat_"+self.settings["sub_model_part_name"].GetString()
        #can not be a child of wall_model_part due to process imposed variables
        self.model_part.CreateSubModelPart(self.contact_model_part_name) 
        self.contact_wall_model_part = self.model_part.GetSubModelPart(self.contact_model_part_name)
        
        print("::[Parametric_Wall]:: -END- ")


    def InitializeSearch(self):

        self.SearchStrategy.InitilizeSearch()        

    def FinalizeSearch(self):

        self.SearchStrategy.FinalizeSearch()
        
    def ExecuteSearch(self):
        
        self.SearchStrategy.ExecuteSearch()
        
