from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import time as timer
import os
import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import main_script as main

sys.path.insert(0,'')
import DEM_explicit_solver_var as DEM_parameters

BaseSolution = main.Solution

class Solution(BaseSolution):

    def __init__(self):
        super(Solution, self).__init__()

    def SetSolverStrategy(self):
        import swimming_sphere_strategy as SolverStrategy
        return SolverStrategy

    def SelectScheme(self):
        scheme = super(Solution, self).SelectScheme()
        if DEM_parameters.IntegrationScheme == 'Hybrid_Bashforth':
            return HybridBashforthScheme()
        elif DEM_parameters.ElementType == "SwimmingNanoParticle":
            return TerminalVelocityScheme()
        else:
            return None

    def SetSolver(self):
        return self.solver_strategy.SwimmingStrategy(self.all_model_parts, self.creator_destructor, self.dem_fem_search, self.scheme, DEM_parameters, self.procedures)

    def ReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        fluid_mp = self.all_model_parts.Get('FluidPart')
        max_node_Id = self.creator_destructor.FindMaxNodeIdInModelPart(fluid_mp)
        max_elem_Id = self.creator_destructor.FindMaxElementIdInModelPart(fluid_mp)
        max_cond_Id = self.creator_destructor.FindMaxConditionIdInModelPart(fluid_mp)
        super(Solution, self).ReadModelParts(max_node_Id + 1, max_elem_Id + 1, max_cond_Id + 1)
