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
    return PostRefiningModeler(main_model_part, meshing_parameters)

class PostRefiningModeler(fluid_mesh_modeler.FluidMeshModeler):
    
    #
    def __init__(self, main_model_part, meshing_parameters): 
        
        self.echo_level        = 1
        self.main_model_part   = main_model_part 
        self.MeshingParameters = meshing_parameters

        print("Construction of the Post Refining Modeler finished")

           
    #
    def InitializeMeshing(self):

        if(self.echo_level>0):
            print("::[fluid_post_refining_modeler]:: -START InitializeMeshing-")

        self.MeshingParameters.InitializeMeshing()

        # set modeler flags: to set options for the mesher (triangle 2D, tetgen 3D)
        # REFINE

        refining_parameters = self.MeshingParameters.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()

        modeler_flags = ""
        modeler_info  = "Refine the domain"

        meshing_options = self.MeshingParameters.GetOptions()

        execution_options = KratosMultiphysics.Flags()

        execution_options.Set(KratosPfem.ModelerUtilities.INITIALIZE_MESHER_INPUT, True)
        execution_options.Set(KratosPfem.ModelerUtilities.FINALIZE_MESHER_INPUT, True)

        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, True)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_ELEMENTS_TO_MESHER, False)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, False)

        if( meshing_options.Is(KratosPfem.ModelerUtilities.CONSTRAINED) ):
            execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, True)

        execution_options.Set(KratosPfem.ModelerUtilities.SELECT_TESSELLATION_ELEMENTS, True)
        execution_options.Set(KratosPfem.ModelerUtilities.KEEP_ISOLATED_NODES, True)


        self.MeshingParameters.SetExecutionOptions(execution_options)

        if( self.dimension == 2 ):
           
            if( refining_options.Is(KratosPfem.ModelerUtilities.REFINE_ADD_NODES) ):
                #"YYJaqrn" "YJq1.4arn" "Jq1.4arn"
                if( meshing_options.Is(KratosPfem.ModelerUtilities.CONSTRAINED) ):
                    modeler_flags = "pYJq1.4arnCQ"  
                else:
                    modeler_flags = "YJq1.4arnQ"

            if( refining_options.Is(KratosPfem.ModelerUtilities.REFINE_INSERT_NODES) ):
                #"riYYJQ" "riYYJQ" "riJQ" "riQ"
                if( meshing_options.Is(KratosPfem.ModelerUtilities.CONSTRAINED) ):
                    modeler_flags = "rinYYJQ"  
                else:
                    #modeler_flags = "nJFMQO4/4"
                    modeler_flags = "nQP"
            
        elif( self.dimension == 3 ):

            if( refining_options.Is(KratosPfem.ModelerUtilities.REFINE_ADD_NODES) ):
                if( meshing_options.Is(KratosPfem.ModelerUtilities.CONSTRAINED) ):
                    modeler_flags = "pMYJq1.4arnCBQF"
                else:
                    modeler_flags = "YJq1.4arnBQF"

            if( refining_options.Is(KratosPfem.ModelerUtilities.REFINE_INSERT_NODES) ):
                if( meshing_options.Is(KratosPfem.ModelerUtilities.CONSTRAINED) ):
                    modeler_flags = "rinYYJBQF"  
                else:
                    # in PSolid "rQYYCCJFu0" 
                    # modeler_flags = "nQMu0"
                    modeler_flags = "nJFu0"
                    #modeler_flags = "QYYCCJF"
                    #modeler_flags = "rQYYCCJF"#PSOLID
                    #modeler_flags = "rinJBQF"

        self.MeshingParameters.SetTessellationFlags(modeler_flags)
        self.MeshingParameters.SetTessellationInfo(modeler_info)


    #
    def SetPreMeshingProcesses(self):

        #remove_mesh_nodes = KratosPfemFluid.RemoveMeshNodesForFluids(self.main_model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPreMeshingProcess(remove_mesh_nodes)
        # no process to start
         #generate_new_nodes  = KratosPfemFluid.GenerateNewNodesBeforeMeshing(self.main_model_part, self.MeshingParameters, self.echo_level)
         #self.mesher.SetPreMeshingProcess(generate_new_nodes)
        pass
     
   #
    def SetPostMeshingProcesses(self):

        # The order set is the order of execution:
        if(self.echo_level>0):
            print("::[fluid_post_refining_modeler]:: -START SetPostMeshingProcesses-")

        #select mesh elements
        #generate_particles  = KratosPfem.GenerateNewNodes(self.main_model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPostMeshingProcess(generate_particles)

        #select mesh elements
        #select_mesh_elements  = KratosPfem.SelectMeshElements(self.main_model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPostMeshingProcess(select_mesh_elements)
        select_mesh_elements  = KratosPfemFluid.SelectMeshElementsForFluids(self.main_model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(select_mesh_elements)


        #rebuild elements
        #rebuild_mesh_elements = KratosPfem.BuildMeshElements(self.main_model_part, self.MeshingParameters, self.echo_level)
        rebuild_mesh_elements = KratosPfem.BuildMeshElements(self.main_model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)

        #rebuild boundary
        rebuild_mesh_boundary = KratosPfem.ReconstructMeshBoundary(self.main_model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_boundary)


    def FinalizeMeshing(self):
        
        print("::[fluid_post_refining_modeler]:: -START FinalizeMeshing-")

        # reset execution flags: to unset the options to be executed in methods and processes
        execution_options = KratosMultiphysics.Flags()
        
        # all flags
        execution_options.Set(KratosPfem.ModelerUtilities.INITIALIZE_MESHER_INPUT, False)
        execution_options.Set(KratosPfem.ModelerUtilities.FINALIZE_MESHER_INPUT, False)

        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, True)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_ELEMENTS_TO_MESHER, False)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, False)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, False)

        execution_options.Set(KratosPfem.ModelerUtilities.SELECT_TESSELLATION_ELEMENTS, True)
        execution_options.Set(KratosPfem.ModelerUtilities.KEEP_ISOLATED_NODES, True)

        self.MeshingParameters.SetExecutionOptions(execution_options)

        self.MeshingParameters.FinalizeMeshing()

    #
