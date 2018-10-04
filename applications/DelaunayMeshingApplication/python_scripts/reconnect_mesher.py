from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mesh mesher (the base class for the mesher derivation)
import mesher

def CreateMesher(main_model_part, meshing_parameters):
    return ReconnectMesher(main_model_part, meshing_parameters)

class ReconnectMesher(mesher.Mesher):

    #
    def __init__(self, main_model_part, meshing_parameters):

        mesher.Mesher.__init__(self, main_model_part, meshing_parameters)

    #
    def InitializeMeshing(self):

        self.MeshingParameters.InitializeMeshing()

        meshing_options = self.MeshingParameters.GetOptions()

        # set execution flags: to set the options to be executed in methods and processes
        execution_options = KratosMultiphysics.Flags()

        execution_options.Set(KratosDelaunay.MesherUtilities.INITIALIZE_MESHER_INPUT, True)
        execution_options.Set(KratosDelaunay.MesherUtilities.FINALIZE_MESHER_INPUT, True)

        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, True)
        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_ELEMENTS_TO_MESHER, False)
        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, False)

        if( meshing_options.Is(KratosDelaunay.MesherUtilities.CONSTRAINED) ):
            execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, True)
            execution_options.Set(KratosDelaunay.MesherUtilities.SELECT_TESSELLATION_ELEMENTS, False)
        else:
            execution_options.Set(KratosDelaunay.MesherUtilities.SELECT_TESSELLATION_ELEMENTS, True)

        execution_options.Set(KratosDelaunay.MesherUtilities.KEEP_ISOLATED_NODES, False)

        self.MeshingParameters.SetExecutionOptions(execution_options)

        # set mesher flags: to set options for the mesher (triangle 2D, tetgen 3D)
        # RECONNECT

        mesher_flags = ""
        mesher_info  = "Reconnect a cloud of points"
        if( self.dimension == 2 ):

            if( meshing_options.Is(KratosDelaunay.MesherUtilities.CONSTRAINED) ):
                mesher_flags = "pnBYYQ"
            else:
                mesher_flags = "nQP"


        elif( self.dimension == 3 ):

            if( meshing_options.Is(KratosDelaunay.MesherUtilities.CONSTRAINED) ):
                mesher_flags = "pnBFMYYQ"
            else:
                mesher_flags = "nBFMQ" # "QJFu0" (1.4.3) "nJFMQO4/4" (1.5.0)

        self.MeshingParameters.SetTessellationFlags(mesher_flags)
        self.MeshingParameters.SetTessellationInfo(mesher_info)

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

    #
    @classmethod
    def _class_prefix(self):
        header = "::[--Reconnect Mesher-]::"
        return header
