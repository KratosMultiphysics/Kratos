from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mesh mesher (the base class for the mesher derivation)
import fluid_mesher

def CreateMesher(main_model_part, meshing_parameters):
    return PreRefiningMesher(main_model_part, meshing_parameters)

class PreRefiningMesher(fluid_mesher.FluidMesher):

    #
    def __init__(self, main_model_part, meshing_parameters):

        self.echo_level        = 1
        self.main_model_part   = main_model_part
        self.MeshingParameters = meshing_parameters

        self.model_part = self.main_model_part
        if( self.main_model_part.Name != self.MeshingParameters.GetSubModelPartName() ):
            self.model_part = self.main_model_part.GetSubModelPart(self.MeshingParameters.GetSubModelPartName())

        print("Construction of the Pre Refining Mesher finished")

    #
    def InitializeMeshing(self):

        if(self.echo_level>0):
            print("::[fluid_pre_refining_mesher]:: -START InitializeMeshing-")

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
                #mesher_flags = "rQYYCCJF"
                #mesher_flags = "nQMu0"
                mesher_flags ="nJQF";
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

        #recover_volume_losses  = KratosPfemFluid.RecoverVolumeLosses(self.model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPreMeshingProcess(recover_volume_losses)

        unactive_peak_elements = False
        unactive_sliver_elements = False
        set_active_flag = KratosPfemFluid.SetActiveFlagMesherProcess(self.main_model_part,unactive_peak_elements,unactive_sliver_elements,self.echo_level)
        self.mesher.SetPreMeshingProcess(set_active_flag)

        inlet_management = KratosPfemFluid.InletManagement(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(inlet_management)

        remove_mesh_nodes = KratosPfemFluid.RemoveMeshNodesForFluids(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(remove_mesh_nodes)

        if( refining_options.Is(KratosDelaunay.MesherUtilities.REFINE_INSERT_NODES) ):
            generate_new_nodes  = KratosPfemFluid.GenerateNewNodesBeforeMeshing(self.model_part, self.MeshingParameters, self.echo_level)
            self.mesher.SetPreMeshingProcess(generate_new_nodes)


    #
    def SetPostMeshingProcesses(self):

        # The order set is the order of execution:
        if(self.echo_level>0):
            print("::[fluid_pre_refining_mesher]:: -START SetPostMeshingProcesses-")


        refining_parameters = self.MeshingParameters.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()


        #select mesh elements
        select_mesh_elements  = KratosPfemFluid.SelectMeshElementsForFluids(self.model_part, self.MeshingParameters, self.echo_level)
        #select_mesh_elements  = KratosDelaunay.SelectMeshElements(self.main_model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(select_mesh_elements)
        #rebuild elements
        #rebuild_mesh_elements = KratosDelaunay.BuildMeshElements(self.main_model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)




        if( refining_options.Is(KratosDelaunay.MesherUtilities.REFINE_ADD_NODES) ):
            select_refine_elements = KratosDelaunay.RefineElementsOnSize(self.model_part, self.MeshingParameters, self.echo_level)
            self.mesher.SetPostMeshingProcess(select_refine_elements)

        #rebuild elements
        rebuild_mesh_elements = KratosDelaunay.GenerateNewElements(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)

        ### rebuild boundary
        ############ choose just one of the following two options: ############        
        ## use this if you want conditions
        ## ATTENTION: this is slow, and must be used together with ModelMeshingWithConditionsForFluids and BuildModelPartBoundary
        #rebuild_mesh_boundary = KratosPfemFluid.GenerateNewConditionsForFluids(self.model_part, self.MeshingParameters, self.echo_level)
        
        ## if you use the following, you will not use/build/compute conditions
        ## ATTENTION: it must be used together with ModelMeshingForFluids and BuildModelPartBoundaryForFluids
        rebuild_mesh_boundary = KratosPfemFluid.BuildMeshBoundaryForFluids(self.model_part, self.MeshingParameters, self.echo_level)
        #######################################################################
        self.mesher.SetPostMeshingProcess(rebuild_mesh_boundary)


        #unactive_peak_elements = False
        #unactive_sliver_elements = False
        #set_active_flag = KratosPfemFluid.SetActiveFlagMesherProcess(self.main_model_part,unactive_peak_elements,unactive_sliver_elements,self.echo_level)
        #self.mesher.SetPostMeshingProcess(set_active_flag)


        #inlet_management = KratosPfemFluid.InletManagement(self.model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPostMeshingProcess(inlet_management)

        #if( refining_options.Is(KratosDelaunay.MesherUtilities.REFINE_INSERT_NODES) ):
            #select_refine_elements = KratosPfemFluid.SetElementsToRefineOnSize(self.model_part, self.MeshingParameters, self.echo_level)
            #self.mesher.SetPostMeshingProcess(select_refine_elements)


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
