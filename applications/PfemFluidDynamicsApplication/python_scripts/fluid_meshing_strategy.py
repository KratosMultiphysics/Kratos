from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

import meshing_strategy

def CreateMeshingStrategy(main_model_part, custom_settings):
    return FluidMeshingStrategy(main_model_part, custom_settings)

class FluidMeshingStrategy(meshing_strategy.MeshingStrategy):

    def SetMeshModelers(self):

        print("::[Fluid Meshing Strategy]:: SET MESH MODELER")

        modelers = []        
        if( self.settings["remesh"].GetBool() and self.settings["refine"].GetBool() ):
            modelers.append("fluid_pre_refining_modeler")
            #modelers.append("fluid_post_refining_modeler")
        elif( self.settings["remesh"].GetBool() ):
            modelers.append("reconnect_modeler")
        elif( self.settings["transfer"].GetBool() ):
            modelers.append("transfer_modeler")
 
        for modeler in modelers:
            meshing_module =__import__(modeler)      
            mesher = meshing_module.CreateMeshModeler(self.main_model_part,self.MeshingParameters) 
            self.mesh_modelers.append(mesher)
  
    #
