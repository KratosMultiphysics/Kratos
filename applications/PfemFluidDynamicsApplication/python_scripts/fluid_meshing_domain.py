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
 
    def ComputeAverageMeshParameters(self):

        self.RefiningParameters.ComputeAndSetMeanVolume(self.main_model_part)
        

    #
    def ComputeInitialAverageMeshParameters(self):

 
        numFluid=0
        mean_nodal_h=0
        for node in self.main_model_part.Nodes:
            if (node.Is(KratosMultiphysics.FLUID)):
                numFluid+=1
                nodal_h=node.GetSolutionStepValue(KratosMultiphysics.NODAL_H)
                mean_nodal_h+=nodal_h 

        mean_nodal_h*=1.0/numFluid;

        print("the mean_nodal_h is  ",mean_nodal_h)
    
        self.RefiningParameters.SetCriticalRadius(mean_nodal_h)
        self.RefiningParameters.SetInitialRadius(mean_nodal_h)
        

    #
