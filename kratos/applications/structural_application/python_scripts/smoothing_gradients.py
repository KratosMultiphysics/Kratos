import os, sys
#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
CheckForPreviousImport()


class Smoothing_Gradients:
    #######################################################################
   
    def __init__(self, model_part, domain_size):

        self.model_part    = model_part
        self.domain_size   = domain_size  
        
    #######################################################################
    def Initialize(self):
        self.smoothing = SmoothingUtility(self.model_part, self.domain_size)   
                 
    #######################################################################   
    def WeightedRecoveryGradients(self):
        self.smoothing.WeightedRecoveryGradients()

     #######################################################################   
    def InterpolatedRecoveryGradients(self):
        self.smoothing.InterpolatedRecoveryGradients()

    def DoubleWeightedRecoveryGradients(self):
         self.smoothing.DoubleWeightedRecoveryGradients() 



    
  

