import KratosMultiphysics as Kratos
from KratosMultiphysics.python_solver import PythonSolver
import KratosMultiphysics.analysis_stage

import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.DEMApplication.DEM_analysis_stage
import KratosMultiphysics.PoromechanicsApplication as Poromechanics
import KratosMultiphysics.PoromechanicsApplication.poromechanics_analysis
import KratosMultiphysics.DemStructuresCouplingApplication as DemStructuresCouplingApplication

class PoroMechanicsCouplingWithDemRadialMultiDofsControlModuleAnalysisStage(Kratos.analysis_stage.AnalysisStage):
    def __init__(self, model, parameters):
        self.parameters = parameters
        self.poromechanics_solution = Poromechanics.poromechanics_analysis.PoromechanicsAnalysis(model, parameters["poromechanics_parameters"])
        self.dem_solution = DEM.DEM_analysis_stage.DEMAnalysisStage(model, parameters["dem_parameters"])
        super().__init__(model, parameters)
        self.effective_stresses_communicator = DemStructuresCouplingApplication.EffectiveStressesCommunicatorUtility(self.poromechanics_solution._GetSolver().main_model_part, self.dem_solution.rigid_face_model_part)
        self.pore_pressure_communicator_utility = DemStructuresCouplingApplication.PorePressureCommunicatorUtility(self.poromechanics_solution._GetSolver().main_model_part, self.dem_solution.spheres_model_part)
        self._CheckCoherentInputs()

    def Initialize(self):
        self.poromechanics_solution.Initialize()
        self.dem_solution.Initialize()
        super().Initialize()
        self.effective_stresses_communicator.Initialize()
        self.pore_pressure_communicator_utility.Initialize()

    def RunSolutionLoop(self):
        while self.poromechanics_solution.KeepAdvancingSolutionLoop():
            self.poromechanics_solution.time = self.poromechanics_solution._GetSolver().AdvanceInTime(self.poromechanics_solution.time)
            self.poromechanics_solution.InitializeSolutionStep()
            self.poromechanics_solution._GetSolver().Predict()
            self.poromechanics_solution._GetSolver().SolveSolutionStep()
            self.poromechanics_solution.FinalizeSolutionStep()
            self.poromechanics_solution.OutputSolutionStep()

            self.effective_stresses_communicator.CopyWallCurrentEffectiveStressesToOldEffectiveStresses()
            self.effective_stresses_communicator.CommunicateCurrentRadialEffectiveStressesToDemWalls()
            self.pore_pressure_communicator_utility.ComputeForceOnParticlesDueToPorePressureGradient()

            time_final_DEM_substepping = self.poromechanics_solution.time
            self.Dt_DEM = self.dem_solution.spheres_model_part.ProcessInfo.GetValue(Kratos.DELTA_TIME)

            #for self.dem_solution.time_dem in self._YieldDEMTime(self.dem_solution.time, time_final_DEM_substepping, self.Dt_DEM):
            #print("Now solving DEM...")
            self.stationarity_checking_is_activated = False
            self.dem_step_counter = 0
            while not self.DEMSolutionIsSteady():
                self.dem_solution.time = self.dem_solution._GetSolver().AdvanceInTime(self.dem_solution.time)
                self.dem_solution.InitializeSolutionStep()
                self.dem_solution._GetSolver().Predict()
                self.dem_solution._GetSolver().SolveSolutionStep()
                self.dem_solution.FinalizeSolutionStep()
                self.dem_solution.OutputSolutionStep()

    def DEMSolutionIsSteady(self):
        tolerance = 1.0e-4
        if self.dem_step_counter == 100:
            self.dem_step_counter = 0

            if not DEM.StationarityChecker().CheckIfVariableIsNullInModelPart(self.dem_solution.spheres_model_part, Kratos.TOTAL_FORCES_X, tolerance):
                self.stationarity_checking_is_activated = True
                return False
            if not DEM.StationarityChecker().CheckIfVariableIsNullInModelPart(self.dem_solution.spheres_model_part, Kratos.TOTAL_FORCES_Y, tolerance):
                self.stationarity_checking_is_activated = True
                return False
            if not DEM.StationarityChecker().CheckIfVariableIsNullInModelPart(self.dem_solution.spheres_model_part, Kratos.TOTAL_FORCES_Z, tolerance):
                self.stationarity_checking_is_activated = True
                return False

            if self.stationarity_checking_is_activated:
                #print("DEM solution is steady...")
                return True
            else:
                return False
        else:
            self.dem_step_counter += 1
            return False

    def InitializeSolutionStep(self):
        self.poromechanics_solution.InitializeSolutionStep()
        super().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.poromechanics_solution.FinalizeSolutionStep()
        super().FinalizeSolutionStep()

    def Finalize(self):
        self.poromechanics_solution.Finalize()
        self.dem_solution.Finalize()
        super().Finalize()

    def _CreateSolver(self):
        return DemPoroMechanicsCouplingSolver(self.poromechanics_solution._GetSolver(), self.dem_solution._GetSolver(), self.parameters)

    def _CheckCoherentInputs(self):
        variables_exported_from_gauss_points_to_nodes = self.parameters["poromechanics_parameters"]["solver_settings"]["gp_to_nodal_variable_list"].GetStringArray()
        if not "WATER_PRESSURE_GRADIENT" in variables_exported_from_gauss_points_to_nodes or not "EFFECTIVE_STRESS_TENSOR" in variables_exported_from_gauss_points_to_nodes:
            raise Exception("Coupling DEM with Poromechanics Error: [\"poromechanics_parameters\"][\"solver_settings\"][\"gp_to_nodal_variable_list\"] must contain \"WATER_PRESSURE_GRADIENT\" and \"EFFECTIVE_STRESS_TENSOR\"\n")

    def _YieldDEMTime(self, current_time, current_time_plus_increment, delta_time):
        current_time += delta_time
        tolerance = 0.0001
        while current_time < (current_time_plus_increment - tolerance * delta_time):
            yield current_time
            current_time += delta_time

        current_time = current_time_plus_increment
        yield current_time

class DemPoroMechanicsCouplingSolver(PythonSolver):
    def __init__(self, poromechanics_solver, dem_solver, parameters):
        self.poromechanics_solver = poromechanics_solver
        self.dem_solver = dem_solver

    def ImportModelPart(self):
        pass

    def GetComputingModelPart(self):
        return self.poromechanics_solver.GetComputingModelPart()

if __name__ == "__main__":
    with open("ProjectParameters.json", 'r') as parameter_file:
        project_parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    PoroMechanicsCouplingWithDemRadialMultiDofsControlModuleAnalysisStage(model, project_parameters).Run()
