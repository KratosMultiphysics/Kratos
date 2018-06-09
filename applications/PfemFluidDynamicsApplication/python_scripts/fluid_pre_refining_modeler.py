from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemApplication as KratosPfem
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mesh modeler (the base class for the modeler derivation)
import fluid_mesh_modeler

def CreateMeshModeler(main_model_part, meshing_parameters):
    return PreRefiningModeler(main_model_part, meshing_parameters)

class PreRefiningModeler(fluid_mesh_modeler.FluidMeshModeler):
    
    #
    def __init__(self, main_model_part, meshing_parameters): 
        
        self.echo_level        = 1
        self.main_model_part   = main_model_part 
        self.MeshingParameters = meshing_parameters
        
        self.model_part = self.main_model_part
        if( self.main_model_part.Name != self.MeshingParameters.GetSubModelPartName() ):
            self.model_part = self.main_model_part.GetSubModelPart(self.MeshingParameters.GetSubModelPartName())
            
        print("Construction of the Pre Refining Modeler finished")
           
    #
    def InitializeMeshing(self):

        if(self.echo_level>0):
            print("::[fluid_pre_refining_modeler]:: -START InitializeMeshing-")

        self.MeshingParameters.InitializeMeshing()
   
        # set execution flags: to set the options to be executed in methods and processes
        meshing_options = self.MeshingParameters.GetOptions()

        execution_options = KratosMultiphysics.Flags()
    
        execution_options.Set(KratosPfem.ModelerUtilities.INITIALIZE_MESHER_INPUT, True)
        execution_options.Set(KratosPfem.ModelerUtilities.FINALIZE_MESHER_INPUT,  True)

        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, True)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_ELEMENTS_TO_MESHER, False)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, False)

        if( meshing_options.Is(KratosPfem.ModelerUtilities.CONSTRAINED) ):
            execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, True)

        execution_options.Set(KratosPfem.ModelerUtilities.SELECT_TESSELLATION_ELEMENTS, True)
        execution_options.Set(KratosPfem.ModelerUtilities.KEEP_ISOLATED_NODES, True)


        self.MeshingParameters.SetExecutionOptions(execution_options)
        
        # set modeler flags: to set options for the mesher (triangle 2D, tetgen 3D)
        # RECONNECT
            
        modeler_flags = ""
        modeler_info  = "Prepare domain for refinement"
        if( self.dimension == 2 ):
           
            if( meshing_options.Is(KratosPfem.ModelerUtilities.CONSTRAINED) ):
                modeler_flags = "pnBYYQ"  
            else:
                modeler_flags = "nQP"

            
        elif( self.dimension == 3 ):

            if( meshing_options.Is(KratosPfem.ModelerUtilities.CONSTRAINED) ):
                modeler_flags = "pnBJFMYYQ"
            else:
                #modeler_flags = "rQYYCCJF"
                #modeler_flags = "nQMu0"
                modeler_flags ="nJFu0";
                #modeler_flags ="VJFu0"; #PSOLID
                #modeler_flags ="rMfjYYaq2.5nQ";
                #modeler_flags = "nJFMQO4/4"

        self.MeshingParameters.SetTessellationFlags(modeler_flags)
        self.MeshingParameters.SetTessellationInfo(modeler_info)


    #
    def SetPreMeshingProcesses(self):

        if(self.echo_level>0):
            print("::[fluid_pre_refining_modeler]:: -START SetPreMeshingProcesses-")

        refining_parameters = self.MeshingParameters.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()

        #recover_volume_losses  = KratosPfemFluid.RecoverVolumeLosses(self.model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPreMeshingProcess(recover_volume_losses)
        
        unactive_peak_elements = False
        unactive_sliver_elements = False
        set_active_flag = KratosPfemFluid.SetActiveFlagProcess(self.main_model_part,unactive_peak_elements,unactive_sliver_elements,self.echo_level)
        self.mesher.SetPreMeshingProcess(set_active_flag)

        inlet_management = KratosPfemFluid.InletManagement(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(inlet_management)

        remove_mesh_nodes = KratosPfemFluid.RemoveMeshNodesForFluids(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(remove_mesh_nodes)

        if( refining_options.Is(KratosPfem.ModelerUtilities.REFINE_INSERT_NODES) ):
            generate_new_nodes  = KratosPfemFluid.GenerateNewNodesBeforeMeshing(self.model_part, self.MeshingParameters, self.echo_level)
            self.mesher.SetPreMeshingProcess(generate_new_nodes)


    #
    def SetPostMeshingProcesses(self):

        # The order set is the order of execution:
        if(self.echo_level>0):
            print("::[fluid_pre_refining_modeler]:: -START SetPostMeshingProcesses-")


        refining_parameters = self.MeshingParameters.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()


        #select mesh elements
        select_mesh_elements  = KratosPfemFluid.SelectMeshElementsForFluids(self.model_part, self.MeshingParameters, self.echo_level)
        #select_mesh_elements  = KratosPfem.SelectMeshElements(self.main_model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(select_mesh_elements)
        #rebuild elements
        #rebuild_mesh_elements = KratosPfem.BuildMeshElements(self.main_model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)




        if( refining_options.Is(KratosPfem.ModelerUtilities.REFINE_ADD_NODES) ):
            select_refine_elements = KratosPfem.SetElementsToRefineOnSize(self.model_part, self.MeshingParameters, self.echo_level)
            self.mesher.SetPostMeshingProcess(select_refine_elements)

      #rebuild elements
        #rebuild_mesh_elements = KratosPfem.BuildMeshElements(self.model_part, self.MeshingParameters, self.echo_level)
        rebuild_mesh_elements = KratosPfem.BuildMeshElements(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)

        #rebuild boundary
        rebuild_mesh_boundary = KratosPfem.BuildMeshBoundary(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_boundary)


        #unactive_peak_elements = False
        #unactive_sliver_elements = False
        #set_active_flag = KratosPfemFluid.SetActiveFlagProcess(self.main_model_part,unactive_peak_elements,unactive_sliver_elements,self.echo_level)
        #self.mesher.SetPostMeshingProcess(set_active_flag)
                

        #inlet_management = KratosPfemFluid.InletManagement(self.model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPostMeshingProcess(inlet_management)

        #if( refining_options.Is(KratosPfem.ModelerUtilities.REFINE_INSERT_NODES) ):
            #select_refine_elements = KratosPfemFluid.SetElementsToRefineOnSize(self.model_part, self.MeshingParameters, self.echo_level)
            #self.mesher.SetPostMeshingProcess(select_refine_elements)


    #
    def FinalizeMeshing(self):
        
        if(self.echo_level>0):
            print("::[fluid_pre_refining_modeler]:: -START FinalizeMeshing-")

        # reset execution flags: to unset the options to be executed in methods and processes
        refining_parameters = self.MeshingParameters.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()

        execution_options = KratosMultiphysics.Flags()

        # set for the post_refining process
        if( refining_options.Is(KratosPfem.ModelerUtilities.REFINE_INSERT_NODES) ):
            execution_options.Set(KratosPfem.ModelerUtilities.INITIALIZE_MESHER_INPUT, False)
            execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, True)
            meshing_options = self.MeshingParameters.GetOptions()

            if( meshing_options.Is(KratosPfem.ModelerUtilities.CONSTRAINED) ):
                execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, True)
                 

        if( refining_options.Is(KratosPfem.ModelerUtilities.REFINE_ADD_NODES) ):
            execution_options.Set(KratosPfem.ModelerUtilities.INITIALIZE_MESHER_INPUT, False)

        execution_options.Set(KratosPfem.ModelerUtilities.FINALIZE_MESHER_INPUT,  True)

        execution_options.Set(KratosPfem.ModelerUtilities.SELECT_TESSELLATION_ELEMENTS, True)
        execution_options.Set(KratosPfem.ModelerUtilities.KEEP_ISOLATED_NODES, True)

        self.MeshingParameters.SetExecutionOptions(execution_options)
