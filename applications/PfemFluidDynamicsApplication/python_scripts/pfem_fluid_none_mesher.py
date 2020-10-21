from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
# Import the mesh mesher (the base class for the mesher derivation)
from KratosMultiphysics.PfemFluidDynamicsApplication import fluid_mesher

def CreateMesher(main_model_part, meshing_parameters):
    return PfemFluidNoneMesher(main_model_part, meshing_parameters)

class PfemFluidNoneMesher(fluid_mesher.FluidMesher):

    #
    def SetPreMeshingProcesses(self):

        if(self.echo_level>0):
            print("::[fluid_mesher]:: -START SetPreMeshingProcesses-")
        
        unactive_peak_elements = False
        unactive_sliver_elements = False
        set_active_flag = KratosPfemFluid.SetActiveFlagMesherProcess(self.main_model_part,unactive_peak_elements,unactive_sliver_elements,self.echo_level)
        self.mesher.SetPreMeshingProcess(set_active_flag)

    #
    def SetPostMeshingProcesses(self):

        # The order set is the order of execution:
        if(self.echo_level>0):
            print("::[fluid_mesher]:: -START SetPostMeshingProcesses-")

        rebuild_mesh_boundary = KratosPfemFluid.BuildMeshBoundaryForFluids(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_boundary)

