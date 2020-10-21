from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

from KratosMultiphysics.DelaunayMeshingApplication import meshing_strategy

from importlib import import_module

def CreateMeshingStrategy(main_model_part, custom_settings):
    return FluidMeshingStrategy(main_model_part, custom_settings)

class FluidMeshingStrategy(meshing_strategy.MeshingStrategy):

    def SetMeshers(self):

        meshers_list = []
        if( self.settings["remesh"].GetBool() and self.settings["refine"].GetBool() ):
            meshers_list.append("KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_complete_mesher")
            #mesher_list.append("fluid_post_refining_mesher")
        elif( self.settings["remesh"].GetBool() and self.settings["transfer"].GetBool() ):
            meshers_list.append("KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_none_mesher")
        elif( self.settings["remesh"].GetBool() ):
            meshers_list.append("KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_keeping_nodes_mesher")

        for mesher in meshers_list:
            full_module_name = mesher
            meshing_module = import_module(full_module_name)
            new_mesher = meshing_module.CreateMesher(self.main_model_part,self.MeshingParameters)
            self.meshers.append(new_mesher)

    #
