from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemApplication as KratosPfem
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

import meshing_strategy

def CreateMeshingStrategy(main_model_part, custom_settings):
    return FluidMeshingStrategy(main_model_part, custom_settings)

class FluidMeshingStrategy(meshing_strategy.MeshingStrategy):

    def GetMeshers(self):

        meshers_list = []
        if( self.settings["remesh"].GetBool() and self.settings["refine"].GetBool() ):
            meshers_list.append("fluid_refining_mesher")
        elif( self.settings["remesh"].GetBool() ):
            meshers_list.append("reconnect_mesher")
        elif( self.settings["transfer"].GetBool() ):
            meshers_list.append("transfer_mesher")

        return meshers_list
    #
