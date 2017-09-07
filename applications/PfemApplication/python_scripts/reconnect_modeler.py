from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemApplication as KratosPfem

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mesh modeler (the base class for the modeler derivation)
import mesh_modeler

def CreateMeshModeler(main_model_part, meshing_parameters):
    return ReconnectModeler(main_model_part, meshing_parameters)

class ReconnectModeler(mesh_modeler.MeshModeler):
    
    #
    def __init__(self, main_model_part, meshing_parameters): 

        mesh_modeler.MeshModeler.__init__(self, main_model_part, meshing_parameters)        

        print("::[Reconnect_Modeler]:: -BUILT-")
  
    #
    def InitializeMeshing(self):

        self.MeshingParameters.InitializeMeshing()

        meshing_options = self.MeshingParameters.GetOptions()
        
        # set execution flags: to set the options to be executed in methods and processes
        execution_options = KratosMultiphysics.Flags()

        execution_options.Set(KratosPfem.ModelerUtilities.INITIALIZE_MESHER_INPUT, True)
        execution_options.Set(KratosPfem.ModelerUtilities.FINALIZE_MESHER_INPUT, True)

        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, True)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_ELEMENTS_TO_MESHER, False)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, False)

        if( meshing_options.Is(KratosPfem.ModelerUtilities.CONSTRAINED) ):
            execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, True)
            execution_options.Set(KratosPfem.ModelerUtilities.SELECT_TESSELLATION_ELEMENTS, False)
        else:
            execution_options.Set(KratosPfem.ModelerUtilities.SELECT_TESSELLATION_ELEMENTS, True)
            
        execution_options.Set(KratosPfem.ModelerUtilities.KEEP_ISOLATED_NODES, False)

        self.MeshingParameters.SetExecutionOptions(execution_options)       
        
        # set modeler flags: to set options for the mesher (triangle 2D, tetgen 3D)
        # RECONNECT
            
        modeler_flags = ""
        modeler_info  = "Reconnect a cloud of points"
        if( self.dimension == 2 ):
           
            if( meshing_options.Is(KratosPfem.ModelerUtilities.CONSTRAINED) ):
                modeler_flags = "pnBYYQ"  
            else:
                modeler_flags = "nQP"

            
        elif( self.dimension == 3 ):

            if( meshing_options.Is(KratosPfem.ModelerUtilities.CONSTRAINED) ):
                modeler_flags = "pnBFMYYQ" 
            else:
                modeler_flags = "nBFMQ" # "QJFu0" (1.4.3) "nJFMQO4/4" (1.5.0)

        self.MeshingParameters.SetTessellationFlags(modeler_flags)
        self.MeshingParameters.SetTessellationInfo(modeler_info)

    #
    def SetPreMeshingProcesses(self):
        #nothing to do: only reconnection    
        pass
        
    #
    def FinalizeMeshing(self):
        
        # reset execution flags: to unset the options to be executed in methods and processes
        execution_options = KratosMultiphysics.Flags() 
        # all false
        self.MeshingParameters.SetExecutionOptions(execution_options)

        self.MeshingParameters.FinalizeMeshing()
