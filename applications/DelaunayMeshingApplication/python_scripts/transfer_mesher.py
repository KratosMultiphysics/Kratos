from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mesh mesher (the base class for the mesher derivation)
import mesher

def CreateMesher(main_model_part, meshing_parameters):
    return TransferMesher(main_model_part, meshing_parameters)

class TransferMesher(mesher.Mesher):

    #
    def __init__(self, main_model_part, meshing_parameters):

        mesher.Mesher.__init__(self, main_model_part, meshing_parameters)

    #
    def InitializeMeshing(self):

        self.MeshingParameters.InitializeMeshing()

        # set execution flags: to set the options to be executed in methods and processes
        execution_options = KratosMultiphysics.Flags()

        execution_options.Set(KratosDelaunay.MesherUtilities.INITIALIZE_MESHER_INPUT, True)
        execution_options.Set(KratosDelaunay.MesherUtilities.FINALIZE_MESHER_INPUT, True)

        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, True)
        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_ELEMENTS_TO_MESHER, True)
        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, True)

        execution_options.Set(KratosDelaunay.MesherUtilities.SELECT_TESSELLATION_ELEMENTS, True)
        execution_options.Set(KratosDelaunay.MesherUtilities.KEEP_ISOLATED_NODES, False)

        self.MeshingParameters.SetExecutionOptions(execution_options)


    #
    def SetPreMeshingProcesses(self):
        #nothing to do: only transfer
        pass


   #
    def SetPostMeshingProcesses(self):

        #nothing to do: only reconnection
        #select mesh elements
        #select_mesh_elements  = KratosDelaunay.SelectElements(self.main_model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPostMeshingProcess(select_mesh_elements)

        #rebuild elements
        rebuild_mesh_elements = KratosDelaunay.GenerateNewElements(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)

        #rebuild boundary
        rebuild_mesh_boundary = KratosDelaunay.GenerateNewConditions(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_boundary)

    #
    @classmethod
    def _class_prefix(self):
        header = "::[--Transfer Mesher--]::"
        return header
