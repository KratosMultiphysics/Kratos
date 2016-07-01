from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateMeshModeler(main_model_part, custom_settings):
    return ReconnectModeler(main_model_part, meshing_parameters, mesh_id)

class ReconnectModeler(mesh_modeler.MeshModeler):
    
    #
    def __init__(self, main_model_part, meshing_parameters, mesh_id): 
        
        self.echo_level        = 0
        self.mesh_id           = mesh_id
        self.main_model_part   = main_model_part 
        self.MeshingParameters = meshing_parameters

        print("Construction of RecoonectModeler finished")
  
    #
    def InitializeMeshing(self):
        
        # set execution flags: to set the options to be executed in methods and processes
        execution_options = KratosMultiphysics.Flags()

        execution_options.Set(ModelerUtilities.SET_NODES, True)
        execution_options.Set(ModelerUtilities.SELECT_ELEMENTS, True)
        execution_options.Set(ModelerUtilities.PASS_ALPHA_SHAPE, True)
        execution_options.Set(ModelerUtilities.DELETE_DATA, True) #delete data at the end

        self.MeshingParameters.SetExecutionOptions(execution_options)
        
        meshing_options = self.MeshingParameters.GetOptions()

        # set modeler flags: to set options for the mesher (triangle 2D, tetgen 3D)
        # RECONNECT
            
        modeler_flags = ""
        modeler_info  = "Reconnect a cloud of points"
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
        
        #nothing to do: only reconnection    


    #
    def FinalizeMeshing(self):
        
        # reset execution flags: to unset the options to be executed in methods and processes
        execution_options = KratosMultiphysics.Flags()
        execution_options.Set(ModelerUtilities.SET_NODES, False)
        execution_options.Set(ModelerUtilities.SET_ELEMENTS, False)
        execution_options.Set(ModelerUtilities.SET_FACES, False)  
        execution_options.Set(ModelerUtilities.SELECT_ELEMENTS, False)
        execution_options.Set(ModelerUtilities.SELECT_NODES, False)
        execution_options.Set(ModelerUtilities.PASS_ALPHA_SHAPE, False)
        execution_options.Set(ModelerUtilities.ENGAGED_NODES, False)

        self.MeshingParameters.SetExecutionOptions(execution_options)
