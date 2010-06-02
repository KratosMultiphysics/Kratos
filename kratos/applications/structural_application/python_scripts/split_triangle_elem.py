import os, sys

#importing the Kratos Library
from Kratos import *
from KratosStructuralApplication import *
#from KratosMeshingApplication import *


class Split_Triangle_Element:
    #######################################################################
   
    def __init__(self, model_part, domain_size):

        self.model_part     =  model_part
        self.domain_size    =  domain_size 
        self.smoothing      =  SmoothingUtility(self.model_part, self.domain_size) 
        self.split_triangle =  IntraFractureTriangle(self.model_part, self.domain_size) 

        self.number_of_avg_elems = 10
        self.number_of_avg_nodes = 10
        self.nodal_neighbour_search = FindNodalNeighboursProcess(self.model_part, self.number_of_avg_elems, self.number_of_avg_nodes)
        
       # self.neighbour_calculator = FindElementalNeighboursProcess(model_part,2,10);
       # self.neighbour_calculator.Execute()

        #self.detect        =  DetectElementsAndNodes(self.model_part, self.domain_size)
        #self.split         =  SplitElements(self.model_part, self.domain_size)
        #self.transfer      =  VariableTransferUtility()
        #self.duplicated    =  False;
        
        
        #self.splitted_nodes_list = []
        #self.map_falilure__list  = []

    #######################################################################
    def Initialize(self):
        self.initialize= False  
    

    def Discrete_Fracture_2D(self):
          
        #self.Smoothing();
        #self.smoothing.SettingNodalValues(self.model_part, self.domain_size)
        #self.transfer.TransferVariablesToNodes(self.model_part, DAMAGE)   
        #self.smoothing.Finalize()  
        #self.detect.DetectNode(self.model_part)  
        for pnode in self.model_part.Nodes:
            pnode.SetValue(SPLIT_NODAL, False)
           



        self.model_part.Nodes[13].SetValue(SPLIT_NODAL, True)
        self.model_part.Nodes[24].SetValue(SPLIT_NODAL, True)
        self.model_part.Nodes[12].SetValue(SPLIT_NODAL, True)

        print self.model_part.Nodes[13].GetValue(SPLIT_NODAL)
        print self.model_part.Nodes[24].GetValue(SPLIT_NODAL)
        #self.pnode = self.model_part.Nodes[13];
        #print self.pnode
        
        map_falilure = Array3()
        map_falilure[0] = 0.00;
        map_falilure[1] = 1.00;
        map_falilure[2] = 0.00;
        
        for pnode in self.model_part.Nodes:  
                 if(pnode != None):   
	            if(pnode.GetValue(SPLIT_NODAL)==True):
                       print pnode 
		       #####detecto el nodo que rompera
		       #self.detect.CalculateMapFailure(pnode,  map_falilure)
		       ###### divido el elemento
                       self.split_triangle.DetectAndSplitElements(self.model_part, pnode, map_falilure)
                       self.nodal_neighbour_search.Execute()
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


