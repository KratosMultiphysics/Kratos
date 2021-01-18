from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

from KratosMultiphysics.DelaunayMeshingApplication import mesher

def CreateMesher(main_model_part, meshing_parameters):
    return FluidMesher(main_model_part, meshing_parameters)

class FluidMesher(mesher.Mesher):

    #
    def __init__(self, main_model_part, meshing_parameters):

        self.echo_level        = 1
        self.main_model_part   = main_model_part
        self.MeshingParameters = meshing_parameters

        self.model_part = self.main_model_part
        if( self.main_model_part.Name != self.MeshingParameters.GetSubModelPartName() ):
            self.model_part = self.main_model_part.GetSubModelPart(self.MeshingParameters.GetSubModelPartName())

        print("Construction of the Refining Mesher finished")

    #
    def InitializeMeshing(self):

        if(self.echo_level>0):
            print("::[fluid_mesher]:: -START InitializeMeshing-")

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
        execution_options.Set(KratosDelaunay.MesherUtilities.REFINE_WALL_CORNER, False)

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
    def FinalizeMeshing(self):
        pass

