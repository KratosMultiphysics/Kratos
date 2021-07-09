import KratosMultiphysics as Kratos
from KratosMultiphysics.python_solver import PythonSolver
import KratosMultiphysics.analysis_stage

import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.PoromechanicsApplication as Poromechanics
import KratosMultiphysics.PoromechanicsApplication.poromechanics_analysis
import KratosMultiphysics.DemStructuresCouplingApplication as DemStructuresCouplingApplication
from sp_2d_rigid_fem_algorithm import DEMAnalysisStage2DSpRigidFem

class PoroMechanicsCouplingWithDemRadialMultiDofsControlModuleAnalysisStage(Kratos.analysis_stage.AnalysisStage):
    def __init__(self, model, parameters):
        self.parameters = parameters
        self.poromechanics_solution = Poromechanics.poromechanics_analysis.PoromechanicsAnalysis(model, parameters["poromechanics_parameters"])
        self.dem_solution = DEMAnalysisStage2DSpRigidFem(model, parameters["dem_parameters"])
        super().__init__(model, parameters)
        self.effective_stresses_communicator = DemStructuresCouplingApplication.EffectiveStressesCommunicatorUtility(self.poromechanics_solution._GetSolver().main_model_part, self.dem_solution.rigid_face_model_part)
        self.pore_pressure_communicator_utility = DemStructuresCouplingApplication.PorePressureCommunicatorUtility(self.poromechanics_solution._GetSolver().main_model_part, self.dem_solution.spheres_model_part)
        self._CheckCoherentInputs()
        self.minimum_number_of_DEM_steps_before_checking_steadiness = 20
        self.stationarity_measuring_tolerance = 1.0e-2

    def Initialize(self):
        self.poromechanics_solution.Initialize()
        self.dem_solution.Initialize()
        super().Initialize()
        self.effective_stresses_communicator.Initialize()
        self.pore_pressure_communicator_utility.Initialize()

    def RunSolutionLoop(self):
        
        while self.poromechanics_solution.KeepAdvancingSolutionLoop():
            print("\n************************************ Solving FEM...\n")
            self.poromechanics_solution.time = self.poromechanics_solution._GetSolver().AdvanceInTime(self.poromechanics_solution.time)
            self.poromechanics_solution.InitializeSolutionStep()
            self.poromechanics_solution._GetSolver().Predict()
            self.poromechanics_solution._GetSolver().SolveSolutionStep()
            self.poromechanics_solution.FinalizeSolutionStep()
            self.poromechanics_solution.OutputSolutionStep()

            self.effective_stresses_communicator.CopyWallCurrentEffectiveStressesToOldEffectiveStresses()
            self.effective_stresses_communicator.CommunicateCurrentRadialEffectiveStressesToDemWalls()
            self.pore_pressure_communicator_utility.ComputeForceOnParticlesDueToPorePressureGradient()

            print("\n************************************ Now solving DEM...\n")
            self.DEM_steps_counter = 0
            self.stationarity_checking_is_activated = False

            while not self.DEMSolutionIsSteady():
                self.dem_solution.time = self.dem_solution._GetSolver().AdvanceInTime(self.dem_solution.time)
                self.dem_solution.InitializeSolutionStep()
                self.dem_solution._GetSolver().Predict()
                self.dem_solution._GetSolver().SolveSolutionStep()
                self.dem_solution.FinalizeSolutionStep()
                self.dem_solution.OutputSolutionStep(self.poromechanics_solution.time)

    def DEMSolutionIsSteady(self):
        
        if self.DEM_steps_counter == self.minimum_number_of_DEM_steps_before_checking_steadiness:

            self.DEM_steps_counter = 0
            
            print("\n************************************ Stationarity will now be checked...\n")
            
            if not DEM.StationarityChecker().CheckIfVariableIsNullInModelPart(self.dem_solution.spheres_model_part, Kratos.TOTAL_FORCES_X, self.stationarity_measuring_tolerance):
                print("\n************************************ F_X is still larger than the tolerance given...\n")
                self.stationarity_checking_is_activated = True
                return False
            if not DEM.StationarityChecker().CheckIfVariableIsNullInModelPart(self.dem_solution.spheres_model_part, Kratos.TOTAL_FORCES_Y, self.stationarity_measuring_tolerance):
                print("\n************************************ F_Y is still larger than the tolerance given...\n")
                self.stationarity_checking_is_activated = True
                return False
            if not DEM.StationarityChecker().CheckIfVariableIsNullInModelPart(self.dem_solution.spheres_model_part, Kratos.TOTAL_FORCES_Z, self.stationarity_measuring_tolerance):
                print("\n************************************ F_Z is still larger than the tolerance given...\n")
                self.stationarity_checking_is_activated = True
                return False

            if self.stationarity_checking_is_activated:
                print("\n************************************ DEM solution is steady!...\n")
                return True
            else:
                return False
        else:
            self.DEM_steps_counter += 1
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
