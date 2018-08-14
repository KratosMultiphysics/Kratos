from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mesher (the base class for the mesher derivation)
import mesher

def CreateMesher(main_model_part, meshing_parameters):
    return PreRefiningMesher(main_model_part, meshing_parameters)

class PreRefiningMesher(mesher.Mesher):

    #
    def __init__(self, main_model_part, meshing_parameters):

        mesher.Mesher.__init__(self, main_model_part, meshing_parameters)

    #
    def InitializeMeshing(self):

        self.MeshingParameters.InitializeMeshing()

        # set execution flags: to set the options to be executed in methods and processes
        meshing_options = self.MeshingParameters.GetOptions()

        execution_options = KratosMultiphysics.Flags()

        execution_options.Set(KratosDelaunay.MesherUtilities.INITIALIZE_MESHER_INPUT, True)
        execution_options.Set(KratosDelaunay.MesherUtilities.FINALIZE_MESHER_INPUT,  False)

        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, True)
        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_ELEMENTS_TO_MESHER, False)
        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, False)

        if( meshing_options.Is(KratosDelaunay.MesherUtilities.CONSTRAINED) ):
            execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, True)

        execution_options.Set(KratosDelaunay.MesherUtilities.SELECT_TESSELLATION_ELEMENTS, True)
        execution_options.Set(KratosDelaunay.MesherUtilities.KEEP_ISOLATED_NODES, False)


        self.MeshingParameters.SetExecutionOptions(execution_options)

        # set mesher flags: to set options for the mesher (triangle 2D, tetgen 3D)
        # RECONNECT

        mesher_flags = ""
        mesher_info  = "Prepare domain for refinement"
        if( self.dimension == 2 ):

            if( meshing_options.Is(KratosDelaunay.MesherUtilities.CONSTRAINED) ):
                mesher_flags = "pnBYYQ"
            else:
                mesher_flags = "nQP"


        elif( self.dimension == 3 ):

            if( meshing_options.Is(KratosDelaunay.MesherUtilities.CONSTRAINED) ):
                mesher_flags = "pnBJFMYYQ"
            else:
                mesher_flags = "nJFMQO4/4"

        self.MeshingParameters.SetTessellationFlags(mesher_flags)
        self.MeshingParameters.SetTessellationInfo(mesher_info)


    #
    def SetPreMeshingProcesses(self):


        # process to refine elements /refine boundary
        refine_mesh_elements = KratosDelaunay.RefineElementsOnThreshold(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(refine_mesh_elements)


        refine_edge_elements = KratosDelaunay.RefineElementsInEdges(self.model_part,self.MeshingParameters,self.echo_level)
        self.mesher.SetPreMeshingProcess(refine_edge_elements)


        refine_mesh_boundary = KratosDelaunay.RefineConditions(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(refine_mesh_boundary)

        # process to remove nodes / remove boundary
        remove_mesh_nodes = KratosDelaunay.RemoveNodes(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(remove_mesh_nodes)


    #
    def SetPostMeshingProcesses(self):

        # The order set is the order of execution:

        refining_parameters = self.MeshingParameters.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()

        #select mesh elements
        select_mesh_elements  = KratosDelaunay.SelectElements(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(select_mesh_elements)

        if( refining_options.Is(KratosDelaunay.MesherUtilities.REFINE_ADD_NODES) ):
            select_refine_elements = KratosDelaunay.RefineElementsOnSize(self.model_part, self.MeshingParameters, self.echo_level)
            self.mesher.SetPostMeshingProcess(select_refine_elements)


        #if( refining_options.Is(KratosDelaunay.MesherUtilities.REFINE_INSERT_NODES) ):
            #insert_nodes = InsertNodes(self.model_part, self.MeshingParameters, self.echo_level)
            #self.mesher.SetPostMeshingProcess(insert_nodes)




    #
    def FinalizeMeshing(self):

        # reset execution flags: to unset the options to be executed in methods and processes
        refining_parameters = self.MeshingParameters.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()

        execution_options = KratosMultiphysics.Flags()

        # set for the post_refining process
        if( refining_options.Is(KratosDelaunay.MesherUtilities.REFINE_INSERT_NODES) ):
            execution_options.Set(KratosDelaunay.MesherUtilities.INITIALIZE_MESHER_INPUT, True)
            execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, True)
            meshing_options = self.MeshingParameters.GetOptions()
            if( meshing_options.Is(KratosDelaunay.MesherUtilities.CONSTRAINED) ):
                execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, True)


        if( refining_options.Is(KratosDelaunay.MesherUtilities.REFINE_ADD_NODES) ):
            execution_options.Set(KratosDelaunay.MesherUtilities.INITIALIZE_MESHER_INPUT, False)

        execution_options.Set(KratosDelaunay.MesherUtilities.FINALIZE_MESHER_INPUT,  True)

        execution_options.Set(KratosDelaunay.MesherUtilities.SELECT_TESSELLATION_ELEMENTS, True)
        execution_options.Set(KratosDelaunay.MesherUtilities.KEEP_ISOLATED_NODES, False)

        self.MeshingParameters.SetExecutionOptions(execution_options)

        self.MeshingParameters.InitializeMeshing() # select tessellation elements is going to be performed again

    #
    @classmethod
    def _class_prefix(self):
        header = "::[----Pre Refining---]::"
        return header
