#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
# Import the mesh mesher (the base class for the mesher derivation)
from KratosMultiphysics.PfemFluidDynamicsApplication import fluid_mesher

def CreateMesher(main_model_part, meshing_parameters):
    return PfemFluidKeepingNodesMesher(main_model_part, meshing_parameters)

class PfemFluidKeepingNodesMesher(fluid_mesher.FluidMesher):

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

    #
    def SetPostMeshingProcesses(self):

        # The order set is the order of execution:
        if(self.echo_level>0):
            print("::[fluid_mesher]:: -START SetPostMeshingProcesses-")

        #select mesh elements
        select_mesh_elements  = KratosPfemFluid.SelectMeshElementsForFluids(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(select_mesh_elements)

        #rebuild elements
        rebuild_mesh_elements = KratosDelaunay.GenerateNewElements(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)

        rebuild_mesh_boundary = KratosPfemFluid.BuildMeshBoundaryForFluids(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_boundary)

