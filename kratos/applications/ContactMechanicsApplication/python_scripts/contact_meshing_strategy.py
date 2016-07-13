from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase

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

        meshing_strategy.MeshingStrategy.Intialize(self,meshing_parameters,domain_size,mesh_id)

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
