from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateMeshModeler(main_model_part, custom_settings):
    return PreRefiningModeler(main_model_part, meshing_parameters, mesh_id)

class PreRefiningModeler(mesh_modeler.MeshModeler):
    
    #
    def __init__(self, main_model_part, meshing_parameters, mesh_id): 
        
        self.echo_level        = 0
        self.mesh_id           = mesh_id
        self.main_model_part   = main_model_part 
        self.MeshingParameters = meshing_parameters

        print("Construction of the Pre Refining Modeler finished")


           
    #
    def InitializeMeshing(self):
        
        # set execution flags: to set the options to be executed in methods and processes
        execution_options = Flags()
        meshing_options = self.MeshingParameters.GetOptions()
    
        self.MeshingParameters.SetNodalIdsFlag(False) #they are not set
        execution_options.Set(ModelerUtilities.SET_NODES, True)

        execution_options.Set(ModelerUtilities.SET_ELEMENTS, False)
        execution_options.Set(ModelerUtilities.SET_NEIGHBOURS, False)

        if( meshing_options.Is(ModelerUtilities.CONSTRAINED) ):
            execution_options.Set(ModelerUtilities.SET_FACES, True)

        execution_options.Set(ModelerUtilities.SELECT_ELEMENTS, True)
        execution_options.Set(ModelerUtilities.SELECT_NODES, False)
        execution_options.Set(ModelerUtilities.PASS_ALPHA_SHAPE, True)
        execution_options.Set(ModelerUtilities.ENGAGED_NODES, False)
        execution_options.Set(ModelerUtilities.DELETE_DATA, False)

        self.MeshingParameters.SetExecutionOptions(execution_options)
        
        # set modeler flags: to set options for the mesher (triangle 2D, tetgen 3D)
        # RECONNECT
            
        modeler_flags = ""
        modeler_info  = "Prepare domain for refinement"
        if( self.domain_size = 2 ):
           
            if( meshing_options.Is(ModelerUtilities.CONSTRAINED) ):
                modeler_flags = "pnBYYQ"  
            elif:
                modeler_flags = "nQP"

            
        elif( self.domain_size = 3 ):

            if( meshing_options.Is(ModelerUtilities.CONSTRAINED) ):
                modeler_flags = "pnBJFMYYQ"
            elif:
                modeler_flags = "nJFMQO4/4"

        self.MeshingParameters.SetTessellationFlags(modeler_flags)
        self.MeshingParameters.SetTessellationInfo(modeler_info)


    #
    def SetPreMeshingProcesses(self):
        

        # process to refine elements /refine boundary
        refine_mesh_elements = SetElementsToRefineOnThreshold(self.main_model_part, self.RefiningParameters, self.mesh_id, self.echo_level)
        self.mesher.SetPreMeshingProcess(refine_mesh_elements)

        #refine_mesh_boundary = RefineMeshBoundary(self.model_part, self.RefiningParameters, self.mesh_id, self.echo_level)            
        #self.mesher.SetPreMeshingProcess(refine_mesh_boundary)

        # process to remove nodes / remove boundary
        remove_mesh_nodes = RemoveMeshNodes(self.main_model_part, self.MeshingParameters,  self.mesh_id, self.echo_level)
        self.mesher.SetPreMeshingProcess(remove_mesh_nodes)
     

    #
    def SetPostMeshingProcesses(self):

        # The order set is the order of execution:

        refining_parameters = self.MeshingVariables.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()

        if( refining_options.Is(ModelerUtilities.ADD_NODES) ):
            select_refine_elements = SetElementsToRefineOnSize(self.model_part, self.MeshingParameters, self.mesh_id, self.echo_level)
            self.mesher.SetPostMeshingProcess(select_refine_elements)


        #if( refining_options.Is(ModelerUtilities.INSERT_NODES) ):
            #insert_nodes = InsertNodes(self.model_part, self.MeshingParameters, self.mesh_id, self.echo_level)
            #self.mesher.SetPostMeshingProcess(insert_nodes)




    #
    def FinalizeMeshing(self):
        
        # reset execution flags: to unset the options to be executed in methods and processes
        execution_options = Flags()

        self.MeshingParameters.SetNodalIdsFlag(True) #they are already set
        execution_options.Set(ModelerUtilities.SET_NODES, False)

        refining_parameters = self.MeshingVariables.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()

        if( refining_options.Is(ModelerUtilities.INSERT_NODES) ):
            self.MeshingParameters.SetNodalIdsFlag(False) #set it again
            execution_options.Set(ModelerUtilities.SET_NODES, True)
        else:
            self.MeshingParameters.SetNodalIdsFlag(True) #they are already set
            execution_options.Set(ModelerUtilities.SET_NODES, False)


        execution_options.Set(ModelerUtilities.SET_ELEMENTS, False)
        execution_options.Set(ModelerUtilities.SET_NEIGHBOURS, False)
        execution_options.Set(ModelerUtilities.SET_FACES, False)  

        # set for the post_refining process
        execution_options.Set(ModelerUtilities.SELECT_ELEMENTS, True)
        execution_options.Set(ModelerUtilities.SELECT_NODES, True)
        execution_options.Set(ModelerUtilities.PASS_ALPHA_SHAPE, True)
        execution_options.Set(ModelerUtilities.ENGAGED_NODES, True) 

        execution_options.Set(ModelerUtilities.DELETE_DATA, True)

        self.MeshingParameters.SetExecutionOptions(execution_options)
