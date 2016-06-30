from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateMeshModeler(main_model_part, custom_settings):
    return PostRefiningModeler(main_model_part, meshing_parameters, mesh_id)

class PostRefiningModeler(mesh_modeler.MeshModeler):
    
    #
    def __init__(self, main_model_part, meshing_parameters, mesh_id): 
        
        self.echo_level        = 0
        self.mesh_id           = mesh_id
        self.main_model_part   = main_model_part 
        self.MeshingParameters = meshing_parameters

        print("Construction of the Post Refining Modeler finished")

           
    #
    def InitializeMeshing(self):
        
        
        # set modeler flags: to set options for the mesher (triangle 2D, tetgen 3D)
        # REFINE

        refining_parameters = self.MeshingVariables.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()

        modeler_flags = ""
        modeler_info  = "Refine the domain"

        if( self.domain_size = 2 ):
           
            if( refining_options.Is(ModelerUtilities.ADD_NODES) ):
                #"YYJaqrn" "YJq1.4arn" "Jq1.4arn"
                if( meshing_options.Is(ModelerUtilities.CONSTRAINED) ):
                    modeler_flags = "pYJq1.4arnCQ"  
                elif:
                    modeler_flags = "YJq1.4arnQ"

            if( refining_options.Is(ModelerUtilities.INSERT_NODES) ):
                #"riYYJQ" "riYYJQ" "riJQ" "riQ"
                if( meshing_options.Is(ModelerUtilities.CONSTRAINED) ):
                    modeler_flags = "rinYYJQ"  
                elif:
                    modeler_flags = "rinJQ"
            
        elif( self.domain_size = 3 ):

            if( refining_options.Is(ModelerUtilities.ADD_NODES) ):
                if( meshing_options.Is(ModelerUtilities.CONSTRAINED) ):
                    modeler_flags = "pMYJq1.4arnCBQF"
                elif:
                    modeler_flags = "YJq1.4arnBQF"

            if( refining_options.Is(ModelerUtilities.INSERT_NODES) ):
                if( meshing_options.Is(ModelerUtilities.CONSTRAINED) ):
                    modeler_flags = "rinYYJBQF"  
                elif:
                    modeler_flags = "rinJBQF"


        self.MeshingParameters.SetTessellationFlags(modeler_flags)
        self.MeshingParameters.SetTessellationInfo(modeler_info)


    #
    def SetPreMeshingProcesses(self):
        
        # no process to start
     

