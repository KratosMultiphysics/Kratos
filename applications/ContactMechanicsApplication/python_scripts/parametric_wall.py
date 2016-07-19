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
    def __init__(self, Model, main_model_part, custom_settings):
        
        self.main_model_part = main_model_part    
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "python_file_name": "parametric_wall",
            "mesh_id": 1,
            "model_part_name" : "Wall Domain",
            "domain_size": 2,
            "echo_level": 1,
            "rigid_body_settings":{
               "rigid_body_element_type": "TranslatoryRigidElement3D1N",
               "fixed_body": true,
               "compute_parameters_from_wall": false,
               "rigid_body_parameters":{
                   "center_of_gravity": [0.0 ,0.0, 0.0],
                   "main_inertias": [0.0, 0.0, 0.0],
                   "principal_axes": [ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ]
               }
            },
            "wall_parameters_list":[{
               "radius": 0.0,
               "center": [0.0, 0.0, 0.0],
               "rake_angle": 0.0,
               "clearance_angle": 0.0,
               "convexity": 0
            }],
            "contact_search_settings":{
               "python_file_name": "contact_search_strategy",
               "search_frequency": 0,
               "contact_parameters":{
                   "contact_condition_type": "PointContactCondition2D1N",
                   "friction_active": false,
                   "friction_law_type": "MorhCoulomb",
                   "variables_of_properties":{
                      "MU_STATIC": 0.3,
                      "MU_DYNAMIC": 0.2,
                      "PENALTY_PARAMETER": 1000
                   }
               }
            }
        }
        """)
        
        ## new node and rigid body element inside the same mesh : boundary conditions also applied
        ## this node and elements must be considered in the computing model part
        ## new contact conditions must be already assembled

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        
        #construct rigid wall // it will contain the array of nodes, array of elements, and the array of conditions
        self.wall_model_part = Model[self.settings["model_part_name"].GetString()]
        self.wall_bounding_box  = KratosContact.SpatialBoundingBox(self.wall_model_part, wall_parameters_list) 

        #construct rigid element // must pass an array of nodes to the element, create a node (CG) and a rigid element set them in the main_model_part, set the node CG as the reference node of the wall_bounding_box, BLOCKED, set in the wall_model_part for imposed movements processes.
        self.rigid_wall_element = KratosContact.CreateRigidBodyElement(self.main_model_part, self.wall_bounding_box, self.settings["rigid_body_settings"])

        #construct the search strategy
        meshing_module = __import__(self.settings["contact_search_strategy"]["python_file_name"].GetString())
        self.SearchStrategy = meshing_module.CreateSearchStrategy(self.main_model_part, self.wall_bounding_box, self.settings["contact_search_settings"])

        print("Construction of Parametric Wall finished")
        

    #### 

    def Initialize(self):

        print("::[Parametric_Wall]:: -START-")
        
        self.domain_size = self.settings["domain_size"].GetInt()
        self.mesh_id     = self.settings["mesh_id"].GetInt()

       
        # Meshing Stratety
        self.SearchStrategy.Initialize()   #in the contact search a contact_mesh or model_part must be created for each parametric wall.

        # create next inside the SearchStrategy: c++
        self.contact_model_part_name =  "contat_"+self.settings["model_part_name"].GetString()
        #can not be a child of wall_model_part due to process imposed variables
        self.main_model_part.CreateSubModelPart(self.contact_model_part_name) 
        self.contact_wall_model_part = self.main_model_part.GetSubModelPart(self.contact_model_part_name)
        
        

        print("::[Parametric_Wall]:: -END- ")


    def InitializeSearch(self):

        self.SearchStrategy.InitilizeSearch()        

    def FinalizeSearch(self):

        self.SearchStrategy.FinalizeSearch()
        
    def ExecuteSearch(self):
        
        self.SearchStrategy.SerchAndBuildContacts()
        
