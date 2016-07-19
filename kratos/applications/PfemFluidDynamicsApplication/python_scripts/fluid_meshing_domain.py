from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

import meshing_domain

def CreateMeshingDomain(main_model_part, custom_settings):
    return FluidMeshingDomain(main_model_part, custom_settings)

class FluidMeshingDomain(meshing_domain.MeshingDomain):
 
    def CheckInitialRadius(self):

        print("::[Mesh Domain]:: CheckInitialRadius -START-")
        
        self.mesh_id     = self.settings["mesh_id"].GetInt()

        self.RefiningParameters.SetInitialRadius(self.main_model_part, self.mesh_id)
        
        print("::[Mesh Domain]:: -END- ")

    #
