from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.ContactMechanicsApplication as KratosContact

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the meshing domain (the base class for the mesher derivation)
import meshing_domain

def CreateMeshingDomain(main_model_part, custom_settings):
    return ContactDomain(main_model_part, custom_settings)

class ContactDomain(meshing_domain.MeshingDomain):

    ##constructor. the constructor shall only take care of storing the settings
    ##and the pointer to the main_model part.
    ##
    ##real construction shall be delayed to the function "Initialize" which
    ##will be called once the mesher is already filled
    def __init__(self, main_model_part, custom_settings):

        self.main_model_part = main_model_part
        self.echo_level      = 1

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
	    "python_module": "contact_domain",
            "model_part_name": "contact_domain",
            "alpha_shape": 1.4,
            "offset_factor": 0.0,
            "meshing_strategy":{
               "python_module": "contact_meshing_strategy",
               "meshing_frequency": 0.0,
               "remesh": true,
               "constrained": true,
               "contact_parameters":{
                   "contact_condition_type": "ContactDomainLM2DCondition",
                   "kratos_module": "KratosMultiphysics.ContactMechanicsApplication",
                   "friction_law_type": "FrictionLaw",
                   "variables_of_properties":{
                       "FRICTION_ACTIVE": false,
                       "MU_STATIC": 0.3,
                       "MU_DYNAMIC": 0.2,
                       "PENALTY_PARAMETER": 1000,
                       "TANGENTIAL_PENALTY_RATIO": 0.1,
                       "TAU_STAB": 1.0
                   }
               }
            },
            "elemental_variables_to_transfer":[ "CAUCHY_STRESS_VECTOR", "DEFORMATION_GRADIENT" ],
            "contact_bodies_list": []
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.RecursivelyValidateAndAssignDefaults(default_settings)

        #construct the solving strategy
        meshing_module = __import__(self.settings["meshing_strategy"]["python_module"].GetString())
        self.MeshingStrategy = meshing_module.CreateMeshingStrategy(self.main_model_part, self.settings["meshing_strategy"])

        self.active_remeshing = False
        if( self.settings["meshing_strategy"]["remesh"].GetBool() ):
            self.active_remeshing = True

        print("::[Contact_Domain]:: (",self.settings["model_part_name"].GetString()," ) -BUILT-")



    ####

    def Initialize(self):

        print("::[Meshing Contact Domain]:: -START-")

        self.dimension = self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION]

        # Set MeshingParameters
        self.SetMeshingParameters()

        # Create contact domain model_part
        if( not self.main_model_part.HasSubModelPart(self.settings["model_part_name"].GetString()) ):
            self.main_model_part.CreateSubModelPart(self.settings["model_part_name"].GetString())

        contact_model_part_names = self.settings["contact_bodies_list"]
        self.contact_parts = KratosContact.StringVector()
        for i in range(contact_model_part_names.size()):
            self.contact_parts.PushBack(contact_model_part_names[i].GetString())

        # transforms list to a std vector
        self.build_contact_model_part = KratosContact.BuildContactModelPart(self.main_model_part, self.MeshingParameters, self.contact_parts, self.echo_level)
        self.build_contact_model_part.Execute()


        # Meshing Stratety
        self.MeshingStrategy.Initialize(self.MeshingParameters, self.dimension)

        print("::[Meshing Contact Domain]:: -END- ")


    ####

    def SetRefiningParameters(self):   #no refine in the contact domain

        # Create RefiningParameters
        self.RefiningParameters = KratosDelaunay.RefiningParameters()
        self.RefiningParameters.Initialize()

        # parameters
        self.RefiningParameters.SetAlphaParameter(self.settings["alpha_shape"].GetDouble())


    def SetMeshingParameters(self):

        # Create MeshingParameters
        meshing_domain.MeshingDomain.SetMeshingParameters(self)


    def ExecuteMeshing(self):

        #if the boundaries has been changed the contact domain has to be updated
        self.build_contact_model_part.Execute()

        self.MeshingStrategy.GenerateMesh()
