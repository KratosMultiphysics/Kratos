from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase
import KratosMultiphysics.ContactMechanicsApplication as KratosContact

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the meshing domain (the base class for the modeler derivation)
import meshing_domain

def CreateMeshingDomain(main_model_part, custom_settings):
    return ContactDomain(main_model_part, custom_settings)

class ContactDomain(mesh_modeler.MeshingDomain):
    
    ##constructor. the constructor shall only take care of storing the settings 
    ##and the pointer to the main_model part.
    ##
    ##real construction shall be delayed to the function "Initialize" which 
    ##will be called once the modeler is already filled
    def __init__(self, main_model_part, custom_settings):
        
        self.main_model_part = main_model_part    
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
	    "python_file_name": "contact_domain",
            "sub_model_part_name": "sub_model_part_name",
            "alpha_shape": 1.4,
            "offset_factor": 0.0,
            "meshing_strategy":{
               "python_file_name": "contact_meshing_strategy",
               "meshing_frequency": 0,
               "remesh": true,
               "constrained": false,
               "contact_parameters":{
                   "contact_condition_type": "ContactDomainLM2DCondition",
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
            },
            "elemental_variables_to_transfer":[ "CAUCHY_STRESS_VECTOR", "DEFORMATION_GRADIENT" ]
        }
        """)
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        #construct the solving strategy
        meshing_module = __import__(self.settings["meshing_strategy"]["python_file_name"].GetString())
        self.MeshingStrategy = meshing_module.CreateMeshingStrategy(self.main_model_part, self.settings["meshing_strategy"])

        print("Construction of Mesh Modeler finished")
        

    #### 

    def Initialize(self):

        print("::[Mesh Contact Domain]:: -START-")
        
        self.domain_size = self.model_part.ProcessInfo[DOMAIN_SIZE]
        self.mesh_id     = 0

        # Set MeshingParameters
        self.SetMeshingParameters()
        
        # Meshing Stratety
        self.MeshingStrategy.Initialize(self.MeshingParameters, self.domain_size, self.mesh_id)
        
        print("::[Mesh Contact Domain]:: -END- ")

        
    ####

    def SetRefiningParameters(self):   #no refine in the contact domain

        # Create RefiningParameters
        self.RefiningParameters = KratosPfemBase.RefiningParameters()
        self.RefiningParameters.Initialize()

        # parameters
        self.RefiningParameters.SetAlphaParameter(self.settings["alpha_shape"].GetDouble())
        

    def SetMeshingParameters(self):
              
        # Create MeshingParameters
        mesh_modeler.MeshingDomain.SetMeshingParameters(self)


    def ExecuteMeshing(self):
        
        self.MeshingStrategy.GenerateMesh()
        
