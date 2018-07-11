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
    return FluidPreRefiningMesher(main_model_part, meshing_parameters)

class FluidPreRefiningMesher(mesher.Mesher):

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
                mesher_flags ="nJFu0";
                #mesher_flags ="VJFu0"; #PSOLID
                #mesher_flags ="rMfjYYaq2.5nQ";
                #mesher_flags = "nJFMQO4/4"

        self.MeshingParameters.SetTessellationFlags(mesher_flags)
        self.MeshingParameters.SetTessellationInfo(mesher_info)

    #
    def SetPreMeshingProcesses(self):

        if(self.echo_level>0):
            print("::[fluid_pre_refining_mesher]:: -START SetPreMeshingProcesses-")

        refining_parameters = self.MeshingParameters.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()

        #recover_volume_losses  = KratosPfem.RecoverVolumeLosses(self.model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPreMeshingProcess(recover_volume_losses)

        unactive_peak_elements = False
        unactive_sliver_elements = False
        set_active_flag = KratosPfem.SetActiveEntities(self.main_model_part,unactive_peak_elements,unactive_sliver_elements,self.echo_level)
        self.mesher.SetPreMeshingProcess(set_active_flag)

        inlet_management = KratosPfem.InletManagement(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(inlet_management)

        remove_mesh_nodes = KratosPfem.RemoveFluidNodes(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(remove_mesh_nodes)

        if( refining_options.Is(KratosDelaunay.MesherUtilities.REFINE_INSERT_NODES) ):
            generate_new_nodes  = KratosPfem.InsertNewNodes(self.model_part, self.MeshingParameters, self.echo_level)
            self.mesher.SetPreMeshingProcess(generate_new_nodes)


    #
    def SetPostMeshingProcesses(self):

        # The order set is the order of execution:
        if(self.echo_level>0):
            print("::[fluid_pre_refining_mesher]:: -START SetPostMeshingProcesses-")

        refining_parameters = self.MeshingParameters.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()

        #select mesh elements
        select_mesh_elements  = KratosPfem.SelectFluidElements(self.model_part, self.MeshingParameters, self.echo_level)
        #select_mesh_elements  = KratosDelaunay.SelectElements(self.main_model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(select_mesh_elements)

        #rebuild elements
        #rebuild_mesh_elements = KratosDelaunay.GenerateNewElements(self.main_model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)

        if( refining_options.Is(KratosDelaunay.MesherUtilities.REFINE_ADD_NODES) ):
            select_refine_elements = KratosDelaunay.RefineElementsOnSize(self.model_part, self.MeshingParameters, self.echo_level)
            self.mesher.SetPostMeshingProcess(select_refine_elements)

        #rebuild elements
        rebuild_mesh_elements = KratosDelaunay.GenerateNewElements(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)

        #rebuild boundary
        rebuild_mesh_boundary = KratosDelaunay.GenerateNewConditions(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_boundary)


    #
    def FinalizeMeshing(self):

        if(self.echo_level>0):
            print("::[fluid_pre_refining_mesher]:: -START FinalizeMeshing-")

        # reset execution flags: to unset the options to be executed in methods and processes
        refining_parameters = self.MeshingParameters.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()

        execution_options = KratosMultiphysics.Flags()

        # set for the post_refining process
        if( refining_options.Is(KratosDelaunay.MesherUtilities.REFINE_INSERT_NODES) ):
            execution_options.Set(KratosDelaunay.MesherUtilities.INITIALIZE_MESHER_INPUT, False)
            execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, True)
            meshing_options = self.MeshingParameters.GetOptions()

            if( meshing_options.Is(KratosDelaunay.MesherUtilities.CONSTRAINED) ):
                execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, True)


        if( refining_options.Is(KratosDelaunay.MesherUtilities.REFINE_ADD_NODES) ):
            execution_options.Set(KratosDelaunay.MesherUtilities.INITIALIZE_MESHER_INPUT, False)

        execution_options.Set(KratosDelaunay.MesherUtilities.FINALIZE_MESHER_INPUT,  True)

        execution_options.Set(KratosDelaunay.MesherUtilities.SELECT_TESSELLATION_ELEMENTS, True)
        execution_options.Set(KratosDelaunay.MesherUtilities.KEEP_ISOLATED_NODES, True)

        self.MeshingParameters.SetExecutionOptions(execution_options)

    #
    def _class_prefix(self):
        header = "::[--Refining Mesher--]::"
        return header
