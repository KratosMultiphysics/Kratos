from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import main_script

BaseAnalysis = main_script.Solution

class Solution(BaseAnalysis):

    def __init__(self, model, pp):
        self.pp = pp
        super(Solution, self).__init__(model)

    def LoadParametersFile(self):
        self.DEM_parameters = self.pp.CFD_DEM

    def GetDefaultInputParameters(self):
        import dem_default_input_parameters
        dem_defaults = dem_default_input_parameters.GetDefaultInputParameters()

        import swimming_dem_default_input_parameters
        only_swimming_defaults = swimming_dem_default_input_parameters.GetDefaultInputParameters()

        for key in only_swimming_defaults.keys():
            dem_defaults.AddValue(key, only_swimming_defaults[key])

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
        translational_scheme = BaseAnalysis.SelectTranslationalScheme(self)
        translational_scheme_name = self.pp.CFD_DEM["TranslationalIntegrationScheme"].GetString()

        if translational_scheme is None:
            if translational_scheme_name == 'Hybrid_Bashforth':
                return HybridBashforthScheme()
            elif translational_scheme_name == "TerminalVelocityScheme":
                return TerminalVelocityScheme()
            else:
                return None
        else:
            return translational_scheme

    def SelectRotationalScheme(self):
        rotational_scheme = BaseAnalysis.SelectRotationalScheme(self)
        translational_scheme_name = self.pp.CFD_DEM["TranslationalIntegrationScheme"].GetString()
        rotational_scheme_name = self.pp.CFD_DEM["RotationalIntegrationScheme"].GetString()

        if rotational_scheme is None:
            if rotational_scheme_name == 'Direct_Integration':
                if translational_scheme_name == 'Hybrid_Bashforth':
                    return HybridBashforthScheme()
                elif translational_scheme_name == 'TerminalVelocityScheme':
                    return TerminalVelocityScheme()
            elif rotational_scheme_name == 'Runge_Kutta':
                return RungeKuttaScheme()
            elif rotational_scheme_name == 'Quaternion_Integration':
                return QuaternionIntegrationScheme()
            else:
                return None
        else:
            return rotational_scheme

    def BaseReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        super(Solution, self).ReadModelParts(max_node_Id, max_elem_Id, max_cond_Id)

    def ReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        self.coupling_analysis.ReadDispersePhaseModelParts()

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
