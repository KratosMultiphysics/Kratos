from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mesh modeler (the base class for the modeler derivation)
import mesh_modeler

def CreateMeshModeler(main_model_part, meshing_parameters):
    return PreRefiningModeler(main_model_part, meshing_parameters)

class PreRefiningModeler(mesh_modeler.MeshModeler):
    
    #
    def __init__(self, main_model_part, meshing_parameters): 
        
        mesh_modeler.MeshModeler.__init__(self, main_model_part, meshing_parameters)        

        print("::[PreRefining_Modeler]:: -BUILT-")
           
    #
    def InitializeMeshing(self):
        
        self.MeshingParameters.InitializeMeshing()

        # set execution flags: to set the options to be executed in methods and processes
        meshing_options = self.MeshingParameters.GetOptions()

        execution_options = KratosMultiphysics.Flags()
    
        execution_options.Set(KratosDelaunay.ModelerUtilities.INITIALIZE_MESHER_INPUT, True)
        execution_options.Set(KratosDelaunay.ModelerUtilities.FINALIZE_MESHER_INPUT,  False)

        execution_options.Set(KratosDelaunay.ModelerUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, True)
        execution_options.Set(KratosDelaunay.ModelerUtilities.TRANSFER_KRATOS_ELEMENTS_TO_MESHER, False)
        execution_options.Set(KratosDelaunay.ModelerUtilities.TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, False)

        if( meshing_options.Is(KratosDelaunay.ModelerUtilities.CONSTRAINED) ):
            execution_options.Set(KratosDelaunay.ModelerUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, True)

        execution_options.Set(KratosDelaunay.ModelerUtilities.SELECT_TESSELLATION_ELEMENTS, True)
        execution_options.Set(KratosDelaunay.ModelerUtilities.KEEP_ISOLATED_NODES, False)


        self.MeshingParameters.SetExecutionOptions(execution_options)
        
        # set modeler flags: to set options for the mesher (triangle 2D, tetgen 3D)
        # RECONNECT
            
        modeler_flags = ""
        modeler_info  = "Prepare domain for refinement"
        if( self.dimension == 2 ):
           
            if( meshing_options.Is(KratosDelaunay.ModelerUtilities.CONSTRAINED) ):
                modeler_flags = "pnBYYQ"  
            else:
                modeler_flags = "nQP"

            
        elif( self.dimension == 3 ):

            if( meshing_options.Is(KratosDelaunay.ModelerUtilities.CONSTRAINED) ):
                modeler_flags = "pnBJFMYYQ"
            else:
                modeler_flags = "nJFMQO4/4"

        self.MeshingParameters.SetTessellationFlags(modeler_flags)
        self.MeshingParameters.SetTessellationInfo(modeler_info)


    #
    def SetPreMeshingProcesses(self):
        

        # process to refine elements /refine boundary
        refine_mesh_elements = KratosDelaunay.SetElementNodesToRefineOnThreshold(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(refine_mesh_elements)


        refine_edge_elements = KratosDelaunay.SetElementEdgesToRefine(self.model_part,self.MeshingParameters,self.echo_level)
        self.mesher.SetPreMeshingProcess(refine_edge_elements)
        

        refine_mesh_boundary = KratosDelaunay.RefineMeshBoundary(self.model_part, self.MeshingParameters, self.echo_level)            
        self.mesher.SetPreMeshingProcess(refine_mesh_boundary)

        # process to remove nodes / remove boundary
        remove_mesh_nodes = KratosDelaunay.RemoveMeshNodes(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(remove_mesh_nodes)
     

    #
    def SetPostMeshingProcesses(self):

        # The order set is the order of execution:

        refining_parameters = self.MeshingParameters.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()

        #select mesh elements
        select_mesh_elements  = KratosDelaunay.SelectMeshElements(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(select_mesh_elements)

        if( refining_options.Is(KratosDelaunay.ModelerUtilities.REFINE_ADD_NODES) ):
            select_refine_elements = KratosDelaunay.SetElementsToRefineOnSize(self.model_part, self.MeshingParameters, self.echo_level)
            self.mesher.SetPostMeshingProcess(select_refine_elements)


        #if( refining_options.Is(KratosDelaunay.ModelerUtilities.REFINE_INSERT_NODES) ):
            #insert_nodes = InsertNodes(self.model_part, self.MeshingParameters, self.echo_level)
            #self.mesher.SetPostMeshingProcess(insert_nodes)




    #
    def FinalizeMeshing(self):
        
        # reset execution flags: to unset the options to be executed in methods and processes
        refining_parameters = self.MeshingParameters.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()

        execution_options = KratosMultiphysics.Flags()

        # set for the post_refining process
        if( refining_options.Is(KratosDelaunay.ModelerUtilities.REFINE_INSERT_NODES) ):
            execution_options.Set(KratosDelaunay.ModelerUtilities.INITIALIZE_MESHER_INPUT, True)
            execution_options.Set(KratosDelaunay.ModelerUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, True)
            meshing_options = self.MeshingParameters.GetOptions()
            if( meshing_options.Is(KratosDelaunay.ModelerUtilities.CONSTRAINED) ):
                execution_options.Set(KratosDelaunay.ModelerUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, True)
                 

        if( refining_options.Is(KratosDelaunay.ModelerUtilities.REFINE_ADD_NODES) ):
            execution_options.Set(KratosDelaunay.ModelerUtilities.INITIALIZE_MESHER_INPUT, False)

        execution_options.Set(KratosDelaunay.ModelerUtilities.FINALIZE_MESHER_INPUT,  True)

        execution_options.Set(KratosDelaunay.ModelerUtilities.SELECT_TESSELLATION_ELEMENTS, True)
        execution_options.Set(KratosDelaunay.ModelerUtilities.KEEP_ISOLATED_NODES, False)

        self.MeshingParameters.SetExecutionOptions(execution_options)
        
        self.MeshingParameters.InitializeMeshing() # select tessellation elements is going to be performed again
