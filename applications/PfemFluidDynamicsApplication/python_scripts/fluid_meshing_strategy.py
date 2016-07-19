from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

import meshing_strategy

def CreateMeshingStrategy(main_model_part, custom_settings):
    return FluidMeshingStrategy(main_model_part, custom_settings)

class FluidMeshingStrategy(meshing_strategy.MeshingStrategy):

    #

    def Initialize(self,meshing_parameters,domain_size,mesh_id):
   
        print("::[fluid_meshing_strategy]:: -START Initialize-")
     
        #parameters
        self.mesh_id = mesh_id

        self.echo_level = 1
        
        #meshing parameters
        self.MeshingParameters = meshing_parameters  
      
        meshing_options = KratosMultiphysics.Flags()
        
        meshing_options.Set(KratosPfemBase.ModelerUtilities.REMESH, self.settings["remesh"].GetBool())
        meshing_options.Set(KratosPfemBase.ModelerUtilities.REFINE, self.settings["refine"].GetBool())
        meshing_options.Set(KratosPfemBase.ModelerUtilities.RECONNECT, self.settings["reconnect"].GetBool())
        meshing_options.Set(KratosPfemBase.ModelerUtilities.TRANSFER, self.settings["transfer"].GetBool())
        meshing_options.Set(KratosPfemBase.ModelerUtilities.CONSTRAINED, self.settings["constrained"].GetBool())
        meshing_options.Set(KratosPfemBase.ModelerUtilities.MESH_SMOOTHING, self.settings["mesh_smoothing"].GetBool())
        meshing_options.Set(KratosPfemBase.ModelerUtilities.VARIABLES_SMOOTHING, self.settings["variables_smoothing"].GetBool())

        self.MeshingParameters.SetOptions(meshing_options)
        self.MeshingParameters.SetReferenceElement(self.settings["reference_element_type"].GetString())
        self.MeshingParameters.SetReferenceCondition(self.settings["reference_condition_type"].GetString())
        
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
        
    #

    def SetMeshModelers(self):

        print("------.......jjjjjjjjjjjjjjjjjjjjjj..------- SetMeshModelers in fluid_meshing_strategy.py")


        modelers = []        
        if( self.settings["remesh"].GetBool() and self.settings["refine"].GetBool() ):
            modelers.append("fluid_pre_refining_modeler")
            modelers.append("fluid_post_refining_modeler")
        elif( self.settings["remesh"].GetBool() ):
            modelers.append("reconnect_modeler")
        elif( self.settings["transfer"].GetBool() ):
            modelers.append("transfer_modeler")
 
        for modeler in modelers:
            meshing_module =__import__(modeler)      
            mesher = meshing_module.CreateMeshModeler(self.main_model_part,self.MeshingParameters,self.mesh_id) 
            self.mesh_modelers.append(mesher)
  
    #

    def InitializeMeshGeneration(self):
        
        self.number_of_elements   = 0
        self.number_of_conditions = 0
        self.number_of_nodes      = 0
        
        info_parameters = self.MeshingParameters.GetInfoParameters()
        info_parameters.Initialize()
        
        self.SetInfo()

        if( self.global_transfer == True ):
            self.MeshDataTransfer.TransferElementalValuesToNodes(self.TransferParameters,self.main_model_part,self.mesh_id)

    #
    def FinalizeMeshGeneration(self):

        self.SetInfo()

        info_parameters    = self.MeshingParameters.GetInfoParameters()
        smoothing_required = info_parameters.CheckMechanicalSmoothing()
        
        refining_parameters = self.MeshingParameters.GetRefiningParameters()
        
        if( self.global_transfer == True ):
            if(smoothing_required):
                #smooth only on selected part based on a threshold variable
                self.MeshDataTransfer.TransferNodalValuesToElementsOnThreshold(self.TransferParameters,refining_parameters,self.main_model_part,self.mesh_id)
            else:
                #smooth all domain
                self.MeshDataTransfer.TransferNodalValuesToElements(self.TransferParameters,self.main_model_part,self.mesh_id)                  
                
        

    #
