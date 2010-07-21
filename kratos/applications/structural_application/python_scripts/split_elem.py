import os, sys
#importing the Kratos Library
from Kratos import *
from KratosStructuralApplication import *
#from KratosMeshingApplication import *


class Split_Elements:
    #######################################################################
   
    def __init__(self, model_part, domain_size):

        self.model_part    =  model_part
        self.domain_size   =  domain_size 
        self.smoothing     =  SmoothingUtility(self.model_part, self.domain_size) 
        self.detect        =  DetectElementsAndNodes(self.model_part, self.domain_size)
        self.split         =  SplitElements(self.model_part, self.domain_size)
        self.transfer      =  VariableTransferUtility()
        self.duplicated    =  False;
        
        
        #self.splitted_nodes_list = []
        #self.map_falilure__list  = []
        self.map_falilure = ZeroVector(3)

    #######################################################################
    def Initialize(self):
        self.initialize= False  
    



    def Find_Elements(self):
          
        #self.Smoothing();
        #self.smoothing.SettingNodalValues(self.model_part, self.domain_size)
        #self.transfer.TransferVariablesToNodes(self.model_part, DAMAGE)   
        #self.smoothing.Finalize()  
        #self.detect.DetectNode(self.model_part)  


        self.pnode = self.model_part.Nodes[13];
        print self.pnode
        
        map_falilure = Array3()
        map_falilure[0] = 0.90;
        map_falilure[1] = 0.30;
        map_falilure[2] = 0.00;
        
        #for pnode in self.model_part.Nodes:
                
                 ###if(pnode != None):   
	         #if (pnode.GetSolutionStepValue(SPLIT_NODAL)==True):
                    ##print pnode 
		    #####detecto el nodo que rompera
		    #self.detect.CalculateMapFailure(pnode,  map_falilure)
		    ###### divido el elemento
        self.detect.DetectElements(self.model_part, self.pnode, map_falilure)
		    ##self.split.Split(self.model_part, pnode, map_falilure)
	  

	##self.detect.Finalize(self.model_part)
	#self.smoothing.SettingNodalValues(self.model_part, self.domain_size)
	#self.smoothing.RecomputeValuesForNewMesh(self.model_part, self.domain_size)
	#self.Smoothing();
	#self.smoothing.Finalize()           
				  
    #################################################################
    def Smoothing(self):
        #self.smoothing.InterpolatedRecoveryGradients(PK2_STRESS_TENSOR, self.model_part, self.domain_size)
        #self.smoothing.InterpolatedRecoveryGradients(GREEN_LAGRANGE_STRAIN_TENSOR, self.model_part, self.domain_size)  
        self.smoothing.WeightedRecoveryGradients(GREEN_LAGRANGE_STRAIN_TENSOR, self.model_part, self.domain_size)
        #self.smoothing.WeightedRecoveryGradients(PK2_STRESS_TENSOR, self.model_part, self.domain_size)
        self.smoothing.DoubleWeightedRecoveryGradients(DAMAGE, self.model_part, self.domain_size)


