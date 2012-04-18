import os, sys
#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
#from KratosMultiphysics.MeshingApplication import *
CheckForPreviousImport()


class Split_Triangle_Elem:
    #######################################################################
   
    def __init__(self, model_part, domain_size):

        self.model_part    =  model_part
        self.domain_size   =  domain_size 
        self.smoothing     =  SmoothingUtility(self.model_part, self.domain_size)
        self.split         =  IntraFractureTriangle(self.model_part, self.domain_size)
                
        
    #######################################################################
    def Initialize(self):
        self.initialize= False  
    

    def Intra_Fracture(self):
        #self.smoothing.SettingNodalValues(self.model_part, self.domain_size)  
        self.smoothing.WeightedRecoveryGradients(DAMAGE,                        NODAL_DAMAGE,  self.model_part, self.domain_size)
        #self.smoothing.WeightedRecoveryGradients(PK2_STRESS_TENSOR,             NODAL_STRESS,  self.model_part, self.domain_size)
        self.smoothing.WeightedRecoveryGradients(GREEN_LAGRANGE_STRAIN_TENSOR,  NODAL_STRAIN , self.model_part, self.domain_size)


        #refine  =  True; 
        #self.split.DetectAndSplitElements(refine) 
	  

	##self.detect.Finalize(self.model_part)
	#self.smoothing.SettingNodalValues(self.model_part, self.domain_size)
	#self.smoothing.RecomputeValuesForNewMesh(self.model_part, self.domain_size)
	#self.Smoothing();
	self.smoothing.Finalize()           
				  
    #################################################################
    def Smoothing(self):
        self.smoothing.SettingNodalValues(self.model_part, self.domain_size)   
        #self.smoothing.InterpolatedRecoveryGradients(PK2_STRESS_TENSOR, self.model_part, self.domain_size)
        #self.smoothing.InterpolatedRecoveryGradients(GREEN_LAGRANGE_STRAIN_TENSOR, self.model_part, self.domain_size)  
        #self.smoothing.WeightedRecoveryGradients(GREEN_LAGRANGE_STRAIN_TENSOR, self.model_part, self.domain_size)
        #self.smoothing.WeightedRecoveryGradients(PK2_STRESS_TENSOR, self.model_part, self.domain_size)
        #self.smoothing.DoubleWeightedRecoveryGradients(DAMAGE, self.model_part, self.domain_size)


