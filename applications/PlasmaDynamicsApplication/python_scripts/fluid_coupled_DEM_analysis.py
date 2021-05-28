from importlib import import_module

import KratosMultiphysics.DEMApplication as DEM
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage as DEMAnalysis

import KratosMultiphysics.SwimmingDEMApplication as SDEM



class FluidCoupledDEMAnalysis(DEMAnalysis):

    def __init__(self, model, project_parameters):       
        self.parameters = project_parameters
        super().__init__(model, project_parameters['dem_parameters'])

        
    def _CreateSolver(self):      
        self.strategy_file_name = self.parameters['dem_parameters']['solver_settings']['strategy'].GetString()
        if self.strategy_file_name == 'swimming_sphere_strategy': 
            SolverStrategy = import_module("KratosMultiphysics.SwimmingDEMApplication" + "." + self.strategy_file_name)
            return SolverStrategy.SwimmingStrategy(self.all_model_parts,
                                                   self.creator_destructor,
                                                   self.dem_fem_search,
                                                   self.parameters,
                                                   self.procedures)    
        else:
            SolverStrategy = import_module("KratosMultiphysics.PlasmaDynamicsApplication" + "." + self.strategy_file_name)
            if self.strategy_file_name == 'plasma_sphere_strategy':       
                return SolverStrategy.PlasmaStrategy(self.all_model_parts,
                                                     self.creator_destructor,
                                                     self.dem_fem_search,
                                                     self.parameters,
                                                     self.procedures)
            elif self.strategy_file_name == 'blood_sphere_strategy':       
                return SolverStrategy.BloodStrategy(self.all_model_parts,
                                                    self.creator_destructor,
                                                    self.dem_fem_search,
                                                    self.parameters,
                                                    self.procedures)

    def SelectTranslationalScheme(self):
        translational_scheme = DEMAnalysis.SelectTranslationalScheme(self)
        translational_scheme_name = self.project_parameters["TranslationalIntegrationScheme"].GetString()

        if translational_scheme is None:

            if translational_scheme_name == 'Forward_Euler':
                return DEM.ForwardEulerScheme()
            elif translational_scheme_name == "Symplectic_Euler":
                return DEM.SymplecticEulerScheme()
            elif translational_scheme_name == "Taylor_Scheme":
                return DEM.TaylorScheme()
            elif translational_scheme_name == "Velocity_Verlet":
                return DEM.VelocityVerletScheme()
            else:
                return None
        else:
            return translational_scheme

    def SelectRotationalScheme(self):
        rotational_scheme = DEMAnalysis.SelectRotationalScheme(self)
        rotational_scheme_name = self.project_parameters["RotationalIntegrationScheme"].GetString()

        if rotational_scheme is None:
            if rotational_scheme_name == 'Direct_Integration':
                return self.SelectTranslationalScheme()
            elif rotational_scheme_name == 'Runge_Kutta':
                return DEM.RungeKuttaScheme()
            elif rotational_scheme_name == 'Quaternion_Integration':
                return DEM.QuaternionIntegrationScheme()
            else:
                return None
        else:
            return rotational_scheme


    def GetParticleHistoryWatcher(self):
        watcher_type = self.parameters["full_particle_history_watcher"].GetString()
        if watcher_type == 'Empty':
            return None
        elif watcher_type == 'ParticlesHistoryWatcher':
            return DEM.ParticlesHistoryWatcher()

    def BaseReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        super().ReadModelParts(max_node_Id, max_elem_Id, max_cond_Id)

    def ReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        self.coupling_analysis.ReadDispersePhaseModelParts()

    def IsTimeToPrintPostProcess(self):
        return self.analytic_data_counter.Tick()

    def SetGraphicalOutput(self):
        pass

    def GraphicalOutputInitialize(self):
        pass

    def PrintResultsForGid(self, time):
        pass

    def GraphicalOutputFinalize(self):
        pass