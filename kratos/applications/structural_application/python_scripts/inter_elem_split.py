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
       
       
        #@ old version
    def Inter_Fracture_Heuristic(self):
        print "INITIALIZING FRACTURE PROCESS "
        is_split = False; 
        is_split = self.split.DetectAndSplitElementsHeuristicFormula(self.model_part)      
	return is_split;   
       
    
    #@ old version
    def Inter_Fracture(self):
        print "INITIALIZING FRACTURE PROCESS "
        self.smoothing.SettingNodalValues(self.model_part, self.domain_size)  
        self.smoothing.WeightedRecoveryGradients(DAMAGE,                        NODAL_DAMAGE,  self.model_part, self.domain_size)
        #self.smoothing.WeightedRecoveryGradients(PK2_STRESS_TENSOR,             NODAL_STRESS,  self.model_part, self.domain_size)
        self.smoothing.WeightedRecoveryGradients(GREEN_LAGRANGE_STRAIN_TENSOR,  NODAL_STRAIN , self.model_part, self.domain_size)
               
        #print "DETECTING AND SPLIT ELEMENTS "
        is_split = False; 
        is_split = self.split.DetectAndSplitElements(self.model_part)   
        
	#self.split.Finalize(self.model_part)
	
	#self.smoothing.RecomputeValuesForNewMesh(self.model_part, self.domain_size)
	#self.Smoothing();
	self.smoothing.Finalize()           
	return is_split;
	
				  
    #################################################################
    #def Smoothing(self):
        #self.smoothing.InterpolatedRecoveryGradients(PK2_STRESS_TENSOR, self.model_part, self.domain_size)
        #self.smoothing.InterpolatedRecoveryGradients(GREEN_LAGRANGE_STRAIN_TENSOR, self.model_part, self.domain_size)  
        #self.smoothing.WeightedRecoveryGradients(GREEN_LAGRANGE_STRAIN_TENSOR, self.model_part, self.domain_size)
        #self.smoothing.WeightedRecoveryGradients(PK2_STRESS_TENSOR, self.model_part, self.domain_size)
        #self.smoothing.DoubleWeightedRecoveryGradients(DAMAGE, self.model_part, self.domain_size)


