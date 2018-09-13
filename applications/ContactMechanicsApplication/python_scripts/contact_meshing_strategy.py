from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.ContactMechanicsApplication as KratosContact

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the meshing strategy (the base class for the mesher derivation)
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
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.RecursivelyValidateAndAssignDefaults(default_settings)

        #print("::[Contact_Mesher_Strategy]:: Construction of Meshing Strategy finished")

    #
    def Initialize(self,meshing_parameters,dimension):

        #parameters
        self.echo_level = 1

        #meshing parameters
        self.MeshingParameters = meshing_parameters

        meshing_options = KratosMultiphysics.Flags()

        meshing_options.Set(KratosDelaunay.MesherUtilities.REMESH, self.settings["remesh"].GetBool())
        meshing_options.Set(KratosDelaunay.MesherUtilities.CONSTRAINED, self.settings["constrained"].GetBool())
        meshing_options.Set(KratosDelaunay.MesherUtilities.CONTACT_SEARCH, True)

        self.MeshingParameters.SetOptions(meshing_options)
        self.MeshingParameters.SetReferenceCondition(self.settings["contact_parameters"]["contact_condition_type"].GetString())

        #set contact properties
        properties = KratosMultiphysics.Properties(0)

        contact_parameters = self.settings["contact_parameters"]

        #build friction law :: pass it as a property parameter
        friction_law_module    = contact_parameters["kratos_module"].GetString()
        friction_law_type_name = contact_parameters["friction_law_type"].GetString()

        #import module if not previously imported
        module = __import__(friction_law_module)
        module_name = (friction_law_module.split("."))[-1]
        FrictionLaw = getattr(getattr(module, module_name), friction_law_type_name)

        friction_law = FrictionLaw()
        #properties have not python interface for this variable type
        #properties.SetValue( KratosContact.FRICTION_LAW, friction_law.Clone() )

        #properties.SetValue(KratosContact.FRICTION_LAW_NAME, friction_law_type_name )
        properties.SetValue(KratosMultiphysics.KratosGlobals.GetVariable("FRICTION_LAW_NAME"), friction_law_type_name )

        contact_variables = contact_parameters["variables_of_properties"]

        #iterators of a json list are not working right now :: must be done by hand:
        properties.SetValue(KratosMultiphysics.KratosGlobals.GetVariable("FRICTION_ACTIVE"), contact_variables["FRICTION_ACTIVE"].GetBool())
        properties.SetValue(KratosMultiphysics.KratosGlobals.GetVariable("MU_STATIC"), contact_variables["MU_STATIC"].GetDouble())
        properties.SetValue(KratosMultiphysics.KratosGlobals.GetVariable("MU_DYNAMIC"), contact_variables["MU_DYNAMIC"].GetDouble())
        properties.SetValue(KratosMultiphysics.KratosGlobals.GetVariable("PENALTY_PARAMETER"), contact_variables["PENALTY_PARAMETER"].GetDouble())
        properties.SetValue(KratosMultiphysics.KratosGlobals.GetVariable("TANGENTIAL_PENALTY_RATIO"), contact_variables["TANGENTIAL_PENALTY_RATIO"].GetDouble())
        properties.SetValue(KratosMultiphysics.KratosGlobals.GetVariable("TAU_STAB"), contact_variables["TAU_STAB"].GetDouble())


        self.MeshingParameters.SetProperties(properties)

        #set variables to global transfer
        self.MeshDataTransfer   = KratosDelaunay.MeshDataTransferUtilities()
        self.TransferParameters = KratosDelaunay.TransferParameters()
        self.global_transfer    = False

        #mesh meshers for the current strategy
        self.meshers = []

        #configure meshers:
        self.SetMeshers();

        for mesher in self.meshers:
            mesher.Initialize(dimension)

        self.number_of_nodes      = 0
        self.number_of_elements   = 0
        self.number_of_conditions = 0


        # prepare model conditions to recieve data
        transfer_parameters = self.MeshingParameters.GetTransferParameters()
        transfer_options = transfer_parameters.GetOptions()

        if( self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False ):
            transfer_options.Set(KratosDelaunay.MeshDataTransferUtilities.INITIALIZE_MASTER_CONDITION, True)
            transfer_parameters.SetOptions(transfer_options)

            self.MeshDataTransfer.TransferBoundaryData(transfer_parameters,self.main_model_part)

        # set flags for the transfer needed for the contact domain
        transfer_options.Set(KratosDelaunay.MeshDataTransferUtilities.INITIALIZE_MASTER_CONDITION, False)
        transfer_options.Set(KratosDelaunay.MeshDataTransferUtilities.MASTER_ELEMENT_TO_MASTER_CONDITION, True)
        transfer_parameters.SetOptions(transfer_options)


    #
    def SetMeshers(self):

        meshers_list = []

        if( self.settings["remesh"].GetBool() ):
            meshers_list.append("contact_mesher")


        for mesher in meshers_list:
            meshing_module =__import__(mesher)
            new_mesher = meshing_module.CreateMesher(self.main_model_part,self.MeshingParameters)
            self.meshers.append(new_mesher)


    #
    def InitializeMeshGeneration(self):

        self.number_of_elements   = 0
        self.number_of_conditions = 0
        self.number_of_nodes      = 0

        info_parameters = self.MeshingParameters.GetInfoParameters()
        info_parameters.Initialize()

        self.SetMeshInfo()

        #Needed to compute effective gaps in the contact domain with lagrangian multipliers
        #in fact not needed anymore, moved to FinalizeSolutionStep of the ContactDomainCondition
        #only recovered when restart

        #if( self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True ):
        #transfer_parameters = self.MeshingParameters.GetTransferParameters()
        #model_part = self.main_model_part.GetSubModelPart(self.MeshingParameters.GetSubModelPartName())
        #self.MeshDataTransfer.TransferBoundaryData(transfer_parameters,model_part)

    #
    def FinalizeMeshGeneration(self):

        self.SetMeshInfo()

    #
    def GenerateMesh(self):

        self.InitializeMeshGeneration()

        for mesher in self.meshers:
            mesher.ExecuteMeshing()

        self.FinalizeMeshGeneration()
