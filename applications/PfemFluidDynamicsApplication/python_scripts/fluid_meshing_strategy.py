from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

import meshing_strategy

def CreateMeshingStrategy(main_model_part, custom_settings):
    return FluidMeshingStrategy(main_model_part, custom_settings)

class FluidMeshingStrategy(meshing_strategy.MeshingStrategy):

    def SetMeshers(self):

        print("::[Fluid Meshing Strategy]:: SET MESHER")

        meshers_list = []
        if( self.settings["remesh"].GetBool() and self.settings["refine"].GetBool() ):
            meshers_list.append("fluid_pre_refining_mesher")
            #mesher_list.append("fluid_post_refining_mesher")
        elif( self.settings["remesh"].GetBool() ):
            meshers_list.append("reconnect_mesher")
        elif( self.settings["transfer"].GetBool() ):
            meshers_list.append("transfer_mesher")

        for mesher in meshers_list:
            meshing_module =__import__(mesher)
            new_mesher = meshing_module.CreateMesher(self.main_model_part,self.MeshingParameters)
            self.meshers.append(new_mesher)

    #
