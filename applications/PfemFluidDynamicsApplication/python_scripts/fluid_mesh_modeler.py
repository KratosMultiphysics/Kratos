from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemApplication as KratosPfem
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

import mesh_modeler

def CreateMeshModeler(main_model_part, meshing_parameters):
    return FluidMeshModeler(main_model_part, meshing_parameters)

class FluidMeshModeler(mesh_modeler.MeshModeler):

    #
    def InitializeMeshing(self):

        self.MeshingParameters.InitializeMeshing()

        
        # set execution flags: to set the options to be executed in methods and processes
        execution_options = KratosMultiphysics.Flags()

        execution_options.Set(KratosPfem.ModelerUtilities.INITIALIZE_MESHER_INPUT, False)
        execution_options.Set(KratosPfem.ModelerUtilities.FINALIZE_MESHER_INPUT, False)

        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, False)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_ELEMENTS_TO_MESHER, False)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, False)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, False)

        execution_options.Set(KratosPfem.ModelerUtilities.SELECT_TESSELLATION_ELEMENTS, False)
        execution_options.Set(KratosPfem.ModelerUtilities.KEEP_ISOLATED_NODES, False)

        self.MeshingParameters.SetExecutionOptions(execution_options)
        
        # set modeler flags: to set options for the mesher (triangle 2D, tetgen 3D)
        if( self.dimension == 2 ):
            pass
            #REFINE
            #ADD NODES
            #to add_nodes automatically and refine the mesh ("q"-quality mesh and "a"-area constraint switches)
            # "YYJaqrn" "YJq1.4arn" "Jq1.4arn"
            #refine
            #modeler_flags = "YJq1.4arnQ" 
            #refine constrained
            #modeler_flags = "pYJq1.4arnCQ"
            
            #INSERT NODES
            #to insert a set of given points and refine the mesh
            # "rinYYJQ" "rinYYJQ" "rinJQ" "rinQ"
            #refine
            #modeler_flags = "rinJQ" 
            #refine constrained
            #modeler_flags = "rinYYJQ"
            
            #refine without adding nodes
            #modeler_flags = "YJrnQ" 
            
            #RECONNECT
            #to reconnect a set of points only
            #modeler_flags = "nQP"
            #constrained
            #modeler_flags = "pnBYYQ"
            
            #BOUNDARY SEARCH
            #to get conectivities, boundaries and neighbours only
            #modeler_flags = "ncEBQ" 
            
        if( self.dimension == 3 ):
            #other flags
            pass
  

    #

    def FinalizeMeshing(self):

        if(self.echo_level>0):
            print("::[fluid_mesh_modeler]:: -END FinalizeMeshing-")
        
        # reset execution flags: to unset the options to be executed in methods and processes
        execution_options = KratosMultiphysics.Flags()
        
        # all flags
        execution_options.Set(KratosPfem.ModelerUtilities.INITIALIZE_MESHER_INPUT, False)
        execution_options.Set(KratosPfem.ModelerUtilities.FINALIZE_MESHER_INPUT, False)

        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, False)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_ELEMENTS_TO_MESHER, False)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, False)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, False)

        execution_options.Set(KratosPfem.ModelerUtilities.SELECT_TESSELLATION_ELEMENTS, False)
        execution_options.Set(KratosPfem.ModelerUtilities.KEEP_ISOLATED_NODES, False)

        self.MeshingParameters.SetExecutionOptions(execution_options)

        self.MeshingParameters.FinalizeMeshing()

    #
