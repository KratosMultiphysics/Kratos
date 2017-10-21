from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import main_script

BaseAlgorithm = main_script.Solution

class Solution(BaseAlgorithm):

    def __init__(self, pp):
        self.pp = pp
        super(Solution,self).__init__()

    def SetSolverStrategy(self):
        import swimming_sphere_strategy as SolverStrategy
        return SolverStrategy

    def SetSolver(self):
        return self.solver_strategy.SwimmingStrategy(self.all_model_parts, self.creator_destructor, self.dem_fem_search, self.scheme, self.pp.CFD_DEM, self.procedures)

    def SelectScheme(self):
        scheme = BaseAlgorithm.SelectScheme(self)
        if scheme == None:
            if self.pp.CFD_DEM["IntegrationScheme"].GetString() == 'Hybrid_Bashforth':
                return HybridBashforthScheme()
            elif self.pp.CFD_DEM["IntegrationScheme"].GetString() == "TerminalVelocityScheme":
                return TerminalVelocityScheme()
            else:
                return None
        else:
            return scheme

    def BaseReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        super(Solution,self).ReadModelParts(max_node_Id, max_elem_Id, max_cond_Id)

    def ReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        self.coupling_algorithm.ReadDispersePhaseModelParts()

    def GetParticleHistoryWatcher(self):
        watcher_type = self.pp.CFD_DEM["full_particle_history_watcher"].GetString()
        
        if watcher_type == 'Empty':
            return None
        elif watcher_type == 'ParticlesHistoryWatcher':
            return ParticlesHistoryWatcher()

    def SetGraphicalOutput(self):
        pass

    def GraphicalOutputInitialize(self):
        pass

    def PrintResultsForGid(self, time):
        pass

    def GraphicalOutputFinalize(self):
        pass
