from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase
import KratosMultiphysics.ContactMechanicsApplication as KratosContact

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the meshing strategy (the base class for the modeler derivation)
import meshing_strategy

def CreateMeshingStrategy(main_model_part, custom_settings):
    return ContactMeshingStrategy(main_model_part, custom_settings)

class ContactMeshingStrategy(meshing_strategy.MeshingStrategy):

    #
    def __init__(self, main_model_part, custom_settings):

        self.main_model_part = main_model_part    
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
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
        }
        """)
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
                
        self.imposed_walls          = []
        self.consider_imposed_walls = False      

        print("Construction of Contact Mesh Modeler finished")
        
    #
    def Initialize(self,meshing_parameters,domain_size,mesh_id):

       #parameters
        self.mesh_id = mesh_id

        self.echo_level = 1
        
        #meshing parameters
        self.MeshingParameters = meshing_parameters  
      
        meshing_options = KratosMultiphysics.Flags()
        
        meshing_options.Set(KratosPfemBase.ModelerUtilities.REMESH, self.settings["remesh"].GetBool())
        meshing_options.Set(KratosPfemBase.ModelerUtilities.CONSTRAINED, self.settings["constrained"].GetBool())

        self.MeshingParameters.SetOptions(meshing_options)
        self.MeshingParameters.SetReferenceCondition(self.settings["contact_paramters"]["contact_condition_type"].GetString())
        
        #set contact properties
        properties = KratosMultiphysics.Properties()
        
        contact_variables = self.settings["contact_parameters"]["variables_of_properties"]

        #build friction law :: pass it as a property parameter
        module = __import__(contact_variables["implemented_in_module"].GetString())      
        friction_law_module    = contact_variables["implemented_in_module"].GetString()
        friction_law_type_name = contact_variables["friction_law_type"].GetString()
        FrictionLaw   = None

        friction_law_type_call = " FrictionLaw = " friction_law_module + friction_law_type_name
        exec(friction_law_type_call)
        
        friction_law = FrictionLaw()

        properties.SetValue(FRICTION_LAW, friction_law.Clone())

        #iterators of a json list are not working right now :: must be done by hand:
        properties.SetValue(KratosMultiphysics.KratosGlobals.GetVariable("FRICTION_ACTIVE"), contact_variables["FRICTION_ACTIVE"])
        properties.SetValue(KratosMultiphysics.KratosGlobals.GetVariable("MU_STATIC"), contact_variables["MU_STATIC"])
        properties.SetValue(KratosMultiphysics.KratosGlobals.GetVariable("MU_DYNAMIC"), contact_variables["MU_DYNAMIC"])
        properties.SetValue(KratosMultiphysics.KratosGlobals.GetVariable("PENALTY_PARAMETER"), contact_variables["PENALTY_PARAMETER"])
        properties.SetValue(KratosMultiphysics.KratosGlobals.GetVariable("TANGENTIAL_PENALTY_RATIO"), contact_variables["TANGENTIAL_PENALTY_RATIO"])
        properties.SetValue(KratosMultiphysics.KratosGlobals.GetVariable("TAU_STAB"), contact_variables["TAU_STAB"])


        self.MeshingParameters.SetProperties(properties)

        #set variables to global transfer
        self.MeshDataTransfer   = KratosPfemBase.MeshDataTransferUtilities()
        self.TransferParameters = KratosPfemBase.TransferParameters()
        self.global_transfer    = False
        if( self.settings["variables_smoothing"].GetBool() == True ):
            self.global_transfer = True
            transfer_variables = self.settings["elemental_variables_to_smooth"]
            #for variable in transfer_variables:
            #    self.TransferParameters.SetVariable( KratosMultiphysics.KratosGlobals.GetVariable( variable.GetString() ) )
            for i in range(0, transfer_variables.size() ):            
                self.TransferParameters.SetVariable(KratosMultiphysics.KratosGlobals.GetVariable(transfer_variables[i].GetString()))
                            

        #mesh modelers for the current strategy
        self.mesh_modelers = []
        
        #configure meshers: 
        self.SetMeshModelers();
        
        for mesher in self.mesh_modelers:
            mesher.Initialize(domain_size)

        self.number_of_nodes      = 0
        self.number_of_elements   = 0
        self.number_of_conditions = 0
        


        # prepare model conditions to recieve data
        transfer_parameters = self.MeshingParameters.GetTransferParameters()
        transfer_options = transfer_parameters.GetOptions()
        
        transfer_options.Set(KratosPfemBase.MeshDataTransferUtilities.INITIALIZE_MASTER_CONDITION, True)
        transfer_parameters.SetOptions(transfer_options)

        self.MeshDataTransfer.TransferBoundaryData(transfer_parameters,self.main_model_part,self.mesh_id)

        # set flags for the transfer needed for the contact domain
        transfer_options.Set(KratosPfemBase.MeshDataTransferUtilities.INITIALIZE_MASTER_CONDITION, False)
        transfer_options.Set(KratosPfemBase.MeshDataTransferUtilities.MASTER_ELEMENT_TO_MASTER_CONDITION, True)
        transfer_parameters.SetOptions(transfer_options)
        

    #
    def SetMeshModelers(self):

        modelers = []       

        if( self.settings["remesh"].GetBool() ):
            modelers.append("contact_modeler")

            
        for modeler in modelers:
            meshing_module =__import__(modeler)      
            mesher = meshing_module.CreateMeshModeler(self.main_model_part,self.MeshingParameters,self.mesh_id) 
            self.mesh_modelers.append(mesher)

        if( self.consider_imposed_walls ):
            for mesher in self.mesh_modelers:
                mesher.SetImposedWalls(self.imposed_walls)

  
    #
    def InitializeMeshGeneration(self):
        
        self.number_of_elements   = 0
        self.number_of_conditions = 0
        self.number_of_nodes      = 0
        
        info_parameters = self.MeshingParameters.GetInfoParameters()
        info_parameters.Initialize()
        
        self.SetInfo()

        transfer_parameters = self.MeshingParameters.GetTransferParameters()
        self.MeshDataTransfer.TransferBoundaryData(transfer_parameters,self.main_model_part,self.mesh_id)

    #
    def FinalizeMeshGeneration(self):

        self.SetInfo()

    #
    def GenerateMesh(self):

        self.InitializeMeshGeneration()

        for mesher in self.mesh_modelers:
            mesher.ExecuteMeshing()

        self.FinalizeMeshGeneration()
