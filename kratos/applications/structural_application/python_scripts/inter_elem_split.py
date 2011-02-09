import os, sys
#importing the Kratos Library
from Kratos import *
from KratosStructuralApplication import *
#from KratosMeshingApplication import *


class Nodal_Split_Elem:
    #######################################################################
   
    def __init__(self, model_part, domain_size):

        self.model_part    =  model_part
        self.domain_size   =  domain_size 
        self.smoothing     =  SmoothingUtility(self.model_part, self.domain_size)
 
        if(self.domain_size==2):  
	    self.split         =  InterFractureTriangle(self.model_part, self.domain_size)
        else:
	    self.split         =  InterFractureTetrahedra(self.model_part, self.domain_size)
        
    #######################################################################
    def Initialize(self):
        self.initialize= False  
       
    

    def Inter_Fracture(self):

        self.smoothing.SettingNodalValues(self.model_part, self.domain_size)  
        #self.smoothing.WeightedRecoveryGradients(DAMAGE, self.model_part, self.domain_size)
        self.smoothing.WeightedRecoveryGradients(PK2_STRESS_TENSOR,            NODAL_STRESS, self.model_part, self.domain_size)
        self.smoothing.WeightedRecoveryGradients(GREEN_LAGRANGE_STRAIN_TENSOR, NODAL_STRAIN ,self.model_part, self.domain_size)
        self.smoothing.Finalize()  
        
        #print "DETECTING AND SPLIT ELEMENTS "
        #self.split.DetectAndSplitElements(self.model_part) 
	  

	##self.detect.Finalize(self.model_part)
	#self.smoothing.SettingNodalValues(self.model_part, self.domain_size)
	#self.smoothing.RecomputeValuesForNewMesh(self.model_part, self.domain_size)
	#self.Smoothing();
	#self.smoothing.Finalize()           
				  
    #################################################################
    #def Smoothing(self):
        #self.smoothing.InterpolatedRecoveryGradients(PK2_STRESS_TENSOR, self.model_part, self.domain_size)
        #self.smoothing.InterpolatedRecoveryGradients(GREEN_LAGRANGE_STRAIN_TENSOR, self.model_part, self.domain_size)  
        #self.smoothing.WeightedRecoveryGradients(GREEN_LAGRANGE_STRAIN_TENSOR, self.model_part, self.domain_size)
        #self.smoothing.WeightedRecoveryGradients(PK2_STRESS_TENSOR, self.model_part, self.domain_size)
        #self.smoothing.DoubleWeightedRecoveryGradients(DAMAGE, self.model_part, self.domain_size)


