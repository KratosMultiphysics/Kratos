from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#Activate it to import in the gdb path:
#import sys
#sys.path.append('/home/cpuigbo/kratos')
#x = input("stopped to allow debug: set breakpoints and press enter to continue");
import time as timer
# Import system python
import os

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication     as KratosSolid
import KratosMultiphysics.ExternalSolversApplication    as KratosSolvers
import KratosMultiphysics.PfemBaseApplication           as KratosPfemBase
import KratosMultiphysics.ContactMechanicsApplication   as KratosContact
import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPfemSolid
import KratosMultiphysics.PfemFluidDynamicsApplication  as KratosPfemFluid
import KratosPfemFluid

class Solution(KratosPfemFluid.Solution):    

        
    def SetGraphicalOutput(self):
        pass      

    def GraphicalOutputExecuteInitialize(self):
        pass
        
    def GraphicalOutputExecuteBeforeSolutionLoop(self):
        pass
                
    def GraphicalOutputExecuteInitializeSolutionStep(self):
        pass
        
    def GraphicalOutputExecuteFinalizeSolutionStep(self):
        pass
        
    def GraphicalOutputPrintOutput(self):
        pass

    def GraphicalOutputExecuteFinalize(self):
        pass                         


if __name__ == "__main__": 
    Solution().Run() 
    
