from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemApplication as KratosPfem
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mesh mesher (the base class for the mesher derivation)
import mesher

def CreateMesher(main_model_part, meshing_parameters):
    return FluidRefiningMesher(main_model_part, meshing_parameters)

class FluidRefiningMesher(mesher.Mesher):

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
        execution_options.Set(KratosDelaunay.MesherUtilities.FINALIZE_MESHER_INPUT,  True)

        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, True)
        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_ELEMENTS_TO_MESHER, False)
        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, False)

        if( meshing_options.Is(KratosDelaunay.MesherUtilities.CONSTRAINED) ):
            execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, True)

        execution_options.Set(KratosDelaunay.MesherUtilities.SELECT_TESSELLATION_ELEMENTS, True)
        execution_options.Set(KratosDelaunay.MesherUtilities.KEEP_ISOLATED_NODES, True)


        self.MeshingParameters.SetExecutionOptions(execution_options)

        # set mesher flags: to set options for the mesher (triangle 2D, tetgen 3D)
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
                #mesher_flags = "rQYYCCJF"
                #mesher_flags = "nQMu0"
                mesher_flags ="nJFQ";
                #mesher_flags ="VJFu0"; #PSOLID
                #mesher_flags ="rMfjYYaq2.5nQ";
                #mesher_flags = "nJFMQO4/4"

        self.MeshingParameters.SetTessellationFlags(mesher_flags)
        self.MeshingParameters.SetTessellationInfo(mesher_info)

    #
    def SetPreMeshingProcesses(self):

        if(self.echo_level>0):
            print(self._class_prefix()+" Set pre meshing processes")

        refining_parameters = self.MeshingParameters.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()

        #insert inlet layer (to be tested)
        #insert_inlet = KratosPfem.InsertInlet(self.model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPreMeshingProcess(insert_inlet)

        #move and remove
        remove_mesh_nodes = KratosPfem.RemoveFluidNodes(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(remove_mesh_nodes)

        if( refining_options.Is(KratosDelaunay.MesherUtilities.REFINE_INSERT_NODES) ):
            insert_fluid_nodes  = KratosPfem.InsertFluidNodes(self.model_part, self.MeshingParameters, self.echo_level)
            self.mesher.SetPreMeshingProcess(insert_fluid_nodes)

        #refine elements that have all nodes in rigid boundary
        refine_edge_elements = KratosPfem.RefineFluidElementsInEdges(self.model_part,self.MeshingParameters,self.echo_level)
        self.mesher.SetPreMeshingProcess(refine_edge_elements)

    #
    def SetPostMeshingProcesses(self):

        # The order set is the order of execution:
        if(self.echo_level>0):
            print(self._class_prefix()+" Set post meshing processes")

        refining_parameters = self.MeshingParameters.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()

        #select mesh elements
        select_mesh_elements  = KratosDelaunay.SelectElements(self.main_model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(select_mesh_elements)

        #rebuild elements
        rebuild_mesh_elements = KratosDelaunay.GenerateNewElements(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)

        #rebuild boundary
        rebuild_mesh_boundary = KratosDelaunay.GenerateNewConditions(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_boundary)


    #
    @classmethod
    def _class_prefix(self):
        header = "::[--Refining Mesher--]::"
        return header
