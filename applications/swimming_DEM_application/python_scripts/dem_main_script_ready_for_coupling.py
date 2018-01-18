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
        super(Solution, self).__init__()
        
    def LoadParametersFile(self):
        self.DEM_parameters = self.pp.CFD_DEM
        
    def GetDefaultInputParameters(self):
        import dem_default_input_parameters
        dem_defaults = dem_default_input_parameters.GetDefaultInputParameters()
    
        import swimming_dem_default_input_parameters
        only_swimming_defaults = swimming_dem_default_input_parameters.GetDefaultInputParameters()
        
        for key in only_swimming_defaults.keys():
            dem_defaults.AddValue(key,only_swimming_defaults[key])
            
        return dem_defaults

    def SetSolverStrategy(self):
        import swimming_sphere_strategy as SolverStrategy
        return SolverStrategy

    def SetSolver(self):
        return self.solver_strategy.SwimmingStrategy(self.all_model_parts,
                                                     self.creator_destructor,
                                                     self.dem_fem_search,
                                                     self.pp.CFD_DEM,
                                                     self.procedures)

    def SelectTranslationalScheme(self):
        translational_scheme = BaseAlgorithm.SelectTranslationalScheme(self)
        if translational_scheme == None:
            if (self.pp.CFD_DEM["TranslationalIntegrationScheme"].GetString() == 'Hybrid_Bashforth'):
                return HybridBashforthScheme()
            elif (self.pp.CFD_DEM["TranslationalIntegrationScheme"].GetString() == "TerminalVelocityScheme"):
                return TerminalVelocityScheme()
            else:
                return None
        else:
            return translational_scheme
        
    def SelectRotationalScheme(self):
        rotational_scheme = BaseAlgorithm.SelectRotationalScheme(self)
        if rotational_scheme == None:
            if (self.pp.CFD_DEM["RotationalIntegrationScheme"].GetString() == 'Direct_Integration'):
                if (self.pp.CFD_DEM["TranslationalIntegrationScheme"].GetString() == 'Hybrid_Bashforth'):
                    return HybridBashforthScheme()
                elif (self.pp.CFD_DEM["TranslationalIntegrationScheme"].GetString() == 'TerminalVelocityScheme'):
                    return TerminalVelocityScheme()
            elif (self.pp.CFD_DEM["RotationalIntegrationScheme"].GetString() == 'Runge_Kutta'):
                return RungeKuttaScheme()
            elif (self.pp.CFD_DEM["RotationalIntegrationScheme"].GetString() == 'Quaternion_Integration'):
                return QuaternionIntegrationScheme()
            else:
                return None
        else:
            return rotational_scheme

    def BaseReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        super(Solution, self).ReadModelParts(max_node_Id, max_elem_Id, max_cond_Id)

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
