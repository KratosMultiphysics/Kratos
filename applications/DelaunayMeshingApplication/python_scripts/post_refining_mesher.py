from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mesher (the base class for the mesher derivation)
import mesher

def CreateMesher(main_model_part, meshing_parameters):
    return PostRefiningMesher(main_model_part, meshing_parameters)

class PostRefiningMesher(mesher.Mesher):

    #
    def __init__(self, main_model_part, meshing_parameters):

        mesher.Mesher.__init__(self, main_model_part, meshing_parameters)

    #
    def InitializeMeshing(self):


        # set mesher flags: to set options for the mesher (triangle 2D, tetgen 3D)
        # REFINE

        refining_parameters = self.MeshingParameters.GetRefiningParameters()
        refining_options = refining_parameters.GetRefiningOptions()

        mesher_flags = ""
        mesher_info  = "Refine the domain"

        meshing_options = self.MeshingParameters.GetOptions()

        if( self.dimension == 2 ):

            if( refining_options.Is(KratosDelaunay.MesherUtilities.REFINE_ADD_NODES) ):
                #"YYJaqrn" "YJq1.4arn" "Jq1.4arn"
                if( meshing_options.Is(KratosDelaunay.MesherUtilities.CONSTRAINED) ):
                    mesher_flags = "pYJq1.4arnCQ"
                else:
                    mesher_flags = "YJq1.4arnQ"

            if( refining_options.Is(KratosDelaunay.MesherUtilities.REFINE_INSERT_NODES) ):
                #"riYYJQ" "riYYJQ" "riJQ" "riQ"
                if( meshing_options.Is(KratosDelaunay.MesherUtilities.CONSTRAINED) ):
                    mesher_flags = "rinYYJQ"
                else:
                    mesher_flags = "rinJQ"

        elif( self.dimension == 3 ):

            if( refining_options.Is(KratosDelaunay.MesherUtilities.REFINE_ADD_NODES) ):
                if( meshing_options.Is(KratosDelaunay.MesherUtilities.CONSTRAINED) ):
                    mesher_flags = "pMYJq1.4arnCBQF"
                else:
                    mesher_flags = "YJq1.4arnBQF"

            if( refining_options.Is(KratosDelaunay.MesherUtilities.REFINE_INSERT_NODES) ):
                if( meshing_options.Is(KratosDelaunay.MesherUtilities.CONSTRAINED) ):
                    mesher_flags = "rinYYJBQF"
                else:
                    mesher_flags = "rinJBQF"

        self.MeshingParameters.SetTessellationFlags(mesher_flags)
        self.MeshingParameters.SetTessellationInfo(mesher_info)


    #
    def SetPreMeshingProcesses(self):
        # no process to start
        pass

    #
    @classmethod
    def _class_prefix(self):
        header = "::[---Post Refining---]::"
        return header
