import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import KratosMultiphysics.DEMApplication.DEM_analysis_stage as DEM_analysis_stage
from importlib import import_module

BaseAnalysis = DEM_analysis_stage.DEMAnalysisStage

class FluidCoupledDEMAnalysisStage(BaseAnalysis):

    def __init__(self, model, project_parameters):
        self.sdem_parameters = project_parameters
        super().__init__(model, project_parameters['dem_parameters'])

    def SetSolverStrategy(self):
        import KratosMultiphysics.SwimmingDEMApplication.swimming_sphere_strategy as SolverStrategy
        return SolverStrategy

    def _CreateSolver(self):
        def SetSolverStrategy():
            strategy_file_name = self.sdem_parameters['dem_parameters']['solver_settings']['strategy'].GetString()
            imported_module = import_module("KratosMultiphysics.SwimmingDEMApplication" + "." + strategy_file_name)
            return imported_module

        return SetSolverStrategy().SwimmingStrategy(self.all_model_parts,
                                                    self.creator_destructor,
                                                    self.dem_fem_search,
                                                    self.sdem_parameters,
                                                    self.procedures)

    def SelectTranslationalScheme(self):
        translational_scheme = BaseAnalysis.SelectTranslationalScheme(self)
        translational_scheme_name = self.project_parameters["TranslationalIntegrationScheme"].GetString()

        # Force terminal velocity scheme
        if translational_scheme is None:
            if translational_scheme_name == 'Hybrid_Bashforth':
                return SDEM.HybridBashforthScheme()
            elif translational_scheme_name == "TerminalVelocityScheme":
                return SDEM.TerminalVelocityScheme()
            else:
                return None
        else:
            return translational_scheme

    def SelectRotationalScheme(self):
        rotational_scheme = BaseAnalysis.SelectRotationalScheme(self)
        translational_scheme_name = self.project_parameters["TranslationalIntegrationScheme"].GetString()
        rotational_scheme_name = self.project_parameters["RotationalIntegrationScheme"].GetString()

        if rotational_scheme is None:
            if rotational_scheme_name == 'Direct_Integration':
                if translational_scheme_name == 'Hybrid_Bashforth':
                    return SDEM.HybridBashforthScheme()
                elif translational_scheme_name == 'TerminalVelocityScheme':
                    return SDEM.TerminalVelocityScheme()
            elif rotational_scheme_name == 'Runge_Kutta':
                return SDEM.RungeKuttaScheme()
            elif rotational_scheme_name == 'Quaternion_Integration':
                return SDEM.QuaternionIntegrationScheme()
            else:
                return None
        else:
            return rotational_scheme

    def BaseReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        super().ReadModelParts(max_node_Id, max_elem_Id, max_cond_Id)

    def ReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        self.coupling_analysis.ReadDispersePhaseModelParts()

    def GetParticleHistoryWatcher(self):
        watcher_type = self.sdem_parameters["full_particle_history_watcher"].GetString()

        if watcher_type == 'Empty':
            return None
        elif watcher_type == 'ParticlesHistoryWatcher':
            return DEM.ParticlesHistoryWatcher()

    def IsTimeToPrintPostProcess(self):
        return self.analytic_data_counter.Tick()

    def PrintResults(self):
        #### GiD IO ##########################################
        if self.IsTimeToPrintPostProcess():
            self.PrintResultsForGid(self.time)
            self.time_old_print = self.time

    # def SetGraphicalOutput(self):
    #     pass

    # def GraphicalOutputInitialize(self):
    #     pass

    def PrintResultsForGid(self, time):
        super().PrintResultsForGid(time)
    #     pass

    # def GraphicalOutputFinalize(self):
    #     pass