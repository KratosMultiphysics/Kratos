from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
# Import the mesh mesher (the base class for the mesher derivation)
from KratosMultiphysics.PfemFluidDynamicsApplication import fluid_mesher

def CreateMesher(main_model_part, meshing_parameters):
    return PfemFluidDoubleMesher(main_model_part, meshing_parameters)

class PfemFluidDoubleMesher(fluid_mesher.FluidMesher):

    #
    def SetPreMeshingProcesses(self):

        if(self.echo_level>0):
            print("::[fluid_mesher]:: -START SetPreMeshingProcesses-")

        #recover_volume_losses  = KratosPfemFluid.RecoverVolumeLosses(self.model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPreMeshingProcess(recover_volume_losses)
        #unactive_peak_elements = False
        #unactive_sliver_elements = False
        #set_active_flag = KratosPfemFluid.SetActiveFlagMesherProcess(self.main_model_part,unactive_peak_elements,unactive_sliver_elements,self.echo_level)
        #self.mesher.SetPreMeshingProcessFirstMesh(set_active_flag)

        inlet_management = KratosPfemFluid.InletManagement(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcessFirstMesh(inlet_management)

        remove_mesh_nodes = KratosPfemFluid.RemoveMeshNodesForFluidsFirstMesh(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcessFirstMesh(remove_mesh_nodes)

        generate_new_nodes  = KratosPfemFluid.GenerateNewNodesBeforeMeshingFirstMesh(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcessFirstMesh(generate_new_nodes)

        #unactive_peak_elements = False
        #unactive_sliver_elements = False
        #set_active_flag = KratosPfemFluid.SetActiveFlagMesherProcess(self.main_model_part,unactive_peak_elements,unactive_sliver_elements,self.echo_level)
        #self.mesher.SetPreMeshingProcessSecondMesh(set_active_flag)

        #inlet_management = KratosPfemFluid.InletManagement(self.model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPreMeshingProcessSecondMesh(inlet_management)

        remove_mesh_nodes = KratosPfemFluid.RemoveMeshNodesForFluidsSecondMesh(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcessSecondMesh(remove_mesh_nodes)

        generate_new_nodes  = KratosPfemFluid.GenerateNewNodesBeforeMeshingSecondMesh(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcessSecondMesh(generate_new_nodes)


    #
    def SetPostMeshingProcesses(self):

        # The order set is the order of execution:
        if(self.echo_level>0):
            print("::[fluid_mesher]:: -START SetPostMeshingProcesses-")

        #select mesh elements
        select_mesh_elements  = KratosPfemFluid.SelectMeshElementsForFluidsFirstMesh(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcessFirstMesh(select_mesh_elements)
        #rebuild elements
        #rebuild_mesh_elements = KratosDelaunay.BuildMeshElements(self.main_model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)

        #rebuild elements
        rebuild_mesh_elements = KratosDelaunay.GenerateNewElements(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcessFirstMesh(rebuild_mesh_elements)

        ### rebuild boundary
        ############ choose just one of the following two options: ############
        ## use this if you want conditions
        ## ATTENTION: this is slow, and must be used together with ModelMeshingWithConditionsForFluids and BuildModelPartBoundary
        #rebuild_mesh_boundary = KratosPfemFluid.GenerateNewConditionsForFluids(self.model_part, self.MeshingParameters, self.echo_level)
        ## if you use the following, you will not use/build/compute conditions
        ## ATTENTION: it must be used together with ModelMeshingForFluids and BuildModelPartBoundaryForFluids
        rebuild_mesh_boundary = KratosPfemFluid.BuildMeshBoundaryForFluids(self.model_part, self.MeshingParameters, self.echo_level)
        #######################################################################
        self.mesher.SetPostMeshingProcessFirstMesh(rebuild_mesh_boundary)

        # The order set is the order of execution:
        if(self.echo_level>0):
            print("::[fluid_mesher]:: -START SetPostMeshingProcesses-")

        #select mesh elements
        select_mesh_elements  = KratosPfemFluid.SelectMeshElementsForFluidsSecondMesh(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcessSecondMesh(select_mesh_elements)
        #rebuild elements
        #rebuild_mesh_elements = KratosDelaunay.BuildMeshElements(self.main_model_part, self.MeshingParameters, self.echo_level)
        #self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)

        #rebuild elements
        rebuild_mesh_elements = KratosDelaunay.GenerateNewElements(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcessSecondMesh(rebuild_mesh_elements)

        ### rebuild boundary
        rebuild_mesh_boundary = KratosPfemFluid.BuildMeshBoundaryForFluids(self.model_part, self.MeshingParameters, self.echo_level)
        #######################################################################
        self.mesher.SetPostMeshingProcessSecondMesh(rebuild_mesh_boundary)

    #
    def ExecuteMeshing(self):

        self.InitializeMeshing()  #set execution flags and mesher flags

        self.mesher.ExecuteMeshingTwice(self.model_part)

        self.FinalizeMeshing()    #set execution flags and mesher flags
