#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
# Import the mesh mesher (the base class for the mesher derivation)
from KratosMultiphysics.PfemFluidDynamicsApplication import fluid_mesher

def CreateMesher(main_model_part, meshing_parameters):
    return CutPfemFluidCompleteMesher(main_model_part, meshing_parameters)

class CutPfemFluidCompleteMesher(fluid_mesher.FluidMesher):

    #
    def SetPreMeshingProcesses(self):

        if(self.echo_level>0):
            print("::[fluid_mesher]:: -START SetPreMeshingProcesses-")

        #recover_volume_losses  = KratosPfemFluid.RecoverVolumeLosses(self.model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPreMeshingProcess(recover_volume_losses)

        unactive_peak_elements = False
        unactive_sliver_elements = False
        set_active_flag = KratosPfemFluid.SetActiveFlagMesherProcess(self.main_model_part,unactive_peak_elements,unactive_sliver_elements,self.echo_level)
        self.mesher.SetPreMeshingProcess(set_active_flag)

        inlet_management = KratosPfemFluid.InletManagement(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(inlet_management)

        remove_mesh_nodes = KratosPfemFluid.RemoveMeshNodesForFluidsCutPfem(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(remove_mesh_nodes)

        generate_new_nodes  = KratosPfemFluid.GenerateNewNodesBeforeMeshingCutPfem(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(generate_new_nodes)

    #
    def SetPostMeshingProcesses(self):

        # The order set is the order of execution:
        if(self.echo_level>0):
            print("::[fluid_mesher]:: -START SetPostMeshingProcesses-")

        #select mesh elements
        select_mesh_elements  = KratosPfemFluid.SelectMeshElementsForFluidsCutPfem(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(select_mesh_elements)
        #rebuild elements
        #rebuild_mesh_elements = KratosDelaunay.BuildMeshElements(self.main_model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)

        #rebuild elements
        rebuild_mesh_elements = KratosDelaunay.GenerateNewElements(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)

        ### rebuild boundary
        ############ choose just one of the following two options: ############
        ## use this if you want conditions
        ## ATTENTION: this is slow, and must be used together with ModelMeshingWithConditionsForFluids and BuildModelPartBoundary
        #rebuild_mesh_boundary = KratosPfemFluid.GenerateNewConditionsForFluids(self.model_part, self.MeshingParameters, self.echo_level)
        ## if you use the following, you will not use/build/compute conditions
        ## ATTENTION: it must be used together with ModelMeshingForFluids and BuildModelPartBoundaryForFluids
        rebuild_mesh_boundary = KratosPfemFluid.BuildMeshBoundaryForFluids(self.model_part, self.MeshingParameters, self.echo_level)
        #######################################################################
        self.mesher.SetPostMeshingProcess(rebuild_mesh_boundary)


