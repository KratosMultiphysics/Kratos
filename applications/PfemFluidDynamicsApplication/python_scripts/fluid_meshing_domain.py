from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemApplication as KratosPfem
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

import meshing_domain

def CreateMeshingDomain(main_model_part, custom_settings):
    return FluidMeshingDomain(main_model_part, custom_settings)

class FluidMeshingDomain(meshing_domain.MeshingDomain):
 
    def ComputeAverageMeshParameters(self):
        
        ModelerUtils = KratosPfem.ModelerUtilities();
        self.domain_volume =  ModelerUtils.ComputeModelPartVolume(self.main_model_part)
        self.element_mean_volume = 0
        
        number_of_elements =  self.main_model_part.NumberOfElements()
        nodes_for_element  =  self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION] + 1

        if(number_of_elements != 0):
            self.element_mean_volume = self.domain_volume/float(number_of_elements*nodes_for_element)

        self.RefiningParameters.SetMeanVolume(self.element_mean_volume)
            
        
    def GetMeanVolume(self):

        return self.element_mean_volume
        
    def GetTotalVolume(self):

        return self.domain_volume

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
        delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        self.main_model_part.ProcessInfo.SetValue(KratosPfemFluid.INITIAL_DELTA_TIME,delta_time)
        self.main_model_part.ProcessInfo.SetValue(KratosPfemFluid.CURRENT_DELTA_TIME,delta_time)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.PREVIOUS_DELTA_TIME,delta_time)
        self.main_model_part.ProcessInfo.SetValue(KratosPfemFluid.TIME_INTERVAL_CHANGED,False)
        

    #
