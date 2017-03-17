from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import time as timer
import os
import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import main_script as main

sys.path.insert(0,'')
import DEM_explicit_solver_var as DEM_parameters

BaseSolution = main.Solution

class Solution(BaseSolution):

    def __init__(self):
        super(Solution,self).__init__()
        
        
        
    def SetSolverStrategy(self):
        import swimming_sphere_strategy as SolverStrategy            
        return SolverStrategy
    
    def SetSolver(self):
        return self.solver_strategy.SwimmingStrategy(self.all_model_parts, self.creator_destructor, self.dem_fem_search, self.scheme, DEM_parameters, self.procedures)
    
    def FinalizeTimeStep(self):
        super(Solution,self).FinalizeTimeStep()

        
    def Finalize():
        super(Solution,self).Finalize()

            
