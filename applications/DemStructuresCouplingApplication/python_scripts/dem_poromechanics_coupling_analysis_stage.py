import KratosMultiphysics as Kratos
from KratosMultiphysics.python_solver import PythonSolver
import KratosMultiphysics.analysis_stage

import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.PoromechanicsApplication as Poromechanics
import KratosMultiphysics.PoromechanicsApplication.poromechanics_analysis
import KratosMultiphysics.DemStructuresCouplingApplication as DemStructuresCouplingApplication
from sp_2d_rigid_fem_algorithm import DEMAnalysisStage2DSpRigidFem
import os
import math

class PoroMechanicsCouplingWithDemRadialMultiDofsControlModuleAnalysisStage(Kratos.analysis_stage.AnalysisStage):
    def __init__(self, model, parameters):
        self.parameters = parameters
        self.poromechanics_solution = Poromechanics.poromechanics_analysis.PoromechanicsAnalysis(model, parameters["poromechanics_parameters"])
        self.dem_solution = DEMAnalysisStage2DSpRigidFem(model, parameters["dem_parameters"])
        super().__init__(model, parameters)
        self.effective_stresses_communicator = DemStructuresCouplingApplication.EffectiveStressesCommunicatorUtility(self.poromechanics_solution._GetSolver().main_model_part, self.dem_solution.rigid_face_model_part)
        self.pore_pressure_communicator_utility = DemStructuresCouplingApplication.PorePressureCommunicatorUtility(self.poromechanics_solution._GetSolver().main_model_part, self.dem_solution.spheres_model_part)
        self._CheckCoherentInputs()

        self.number_of_DEM_steps_between_steadiness_checks = 25
        self.stationarity_measuring_tolerance = 5e-3

        file1 = os.path.join(os.getcwd(), "poro_solution_time_vs_sp_chunks.txt")
        file2 = os.path.join(os.getcwd(), "poro_solution_time_vs_sp_standard.txt")
        file3 = os.path.join(os.getcwd(), "poro_solution_time_vs_dem_time_vs_dems_steps.txt")
        file4 = os.path.join(os.getcwd(), "radial_normal_and_smoothed_reaction_stresses_values.txt")
        if os.path.exists(file1):
            os.remove(file1)
        if os.path.exists(file2):
            os.remove(file2)
        if os.path.exists(file3):
            os.remove(file3)
        if os.path.exists(file4):
            os.remove(file4)

    def Initialize(self):
        self.poromechanics_solution.Initialize()
        self.dem_solution.Initialize()
        super().Initialize()
        self.effective_stresses_communicator.Initialize()
        self.pore_pressure_communicator_utility.Initialize()

    def RunSolutionLoop(self):
        self.stationarity_checking_is_activated = False

        while self.poromechanics_solution.KeepAdvancingSolutionLoop():
            print("\n************************************ Solving FEM...\n", flush=True)
            self.poromechanics_solution.time = self.poromechanics_solution._GetSolver().AdvanceInTime(self.poromechanics_solution.time)
            self.poromechanics_solution.InitializeSolutionStep()
            self.poromechanics_solution._GetSolver().Predict()
            self.poromechanics_solution._GetSolver().SolveSolutionStep()
            self.poromechanics_solution.FinalizeSolutionStep()
            self.poromechanics_solution.OutputSolutionStep()

            self.effective_stresses_communicator.CopyWallCurrentEffectiveStressesToOldEffectiveStresses()
            self.effective_stresses_communicator.CommunicateCurrentRadialEffectiveStressesToDemWalls()
            self.pore_pressure_communicator_utility.ComputeForceOnParticlesDueToPorePressureGradient()

            print("\n************************************ Now solving DEM...\n", flush=True)
            self.DEM_steps_counter = 0
            self.dem_time_steps_per_fem_time_step = 0
            self.dem_solution.time = self.poromechanics_solution.time

            while not self.DEMSolutionIsSteady():
                self.dem_solution.time = self.dem_solution._GetSolver().AdvanceInTime(self.dem_solution.time)
                self.dem_solution.InitializeSolutionStep()
                self.dem_solution._GetSolver().Predict()
                self.dem_solution._GetSolver().SolveSolutionStep()
                self.dem_solution.FinalizeSolutionStep()
                self.AdditionalDEMOperationsToFinalizeDEMSolutionStep()
                self.dem_solution.OutputSolutionStep()

    def AdditionalDEMOperationsToFinalizeDEMSolutionStep(self):

        if self.dem_solution.IsTimeToPrintPostProcess():

            filename = 'radial_normal_and_smoothed_reaction_stresses_values_' + str(self.dem_solution.time) + '.txt'

            try:
                os.remove(filename)
            except OSError:
                pass

            radial_normal_and_smoothed_reaction_stresses_file = open(filename, 'a')

            average_ratio = 0.0

            for node in self.dem_solution.rigid_face_model_part.Nodes:
                total_radial_normal_stress_value = node.GetValue(DEM.RADIAL_NORMAL_STRESS_COMPONENT)
                A =  node.GetValue(DEM.SMOOTHED_REACTION_STRESS_X)
                B =  node.GetValue(DEM.SMOOTHED_REACTION_STRESS_Y)
                C =  node.GetValue(DEM.SMOOTHED_REACTION_STRESS_Z)
                total_smoothed_reaction_stress_value = math.sqrt(A*A + B*B + C*C)
                if abs(total_radial_normal_stress_value) > 0.0:
                    this_node_applied_stress_ratio = -total_smoothed_reaction_stress_value / total_radial_normal_stress_value
                else:
                    this_node_applied_stress_ratio = 0.0
                average_ratio += this_node_applied_stress_ratio
                displ = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
                displ_norm = math.sqrt(displ[0]*displ[0]+displ[1]*displ[1]+displ[2]*displ[2])
                vel = node.GetSolutionStepValue(Kratos.VELOCITY)
                vel_norm = math.sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2])
                radial_normal_and_smoothed_reaction_stresses_file.write(str(node.Id) + " " + str(node.X) + " " + str(node.Y) + " " + str(node.Z) + " " + str(total_radial_normal_stress_value) + " " + str(total_smoothed_reaction_stress_value) + " " + str(displ_norm) + " " + str(vel_norm) + '\n')

            radial_normal_and_smoothed_reaction_stresses_file.close()
            average_ratio = average_ratio / self.dem_solution.rigid_face_model_part.NumberOfNodes()
            self.dem_solution.KratosPrintInfo("-------- Averaged loading ratio (average quotient current stress/ target stress) is  " + str(average_ratio) + " ------------ \n")

            total_average_target_radial_normal_stress = 0.0
            total_forces_modulus = 0.0

            for node in self.dem_solution.rigid_face_model_part.Nodes:

                total_radial_normal_stress_value = node.GetValue(DEM.RADIAL_NORMAL_STRESS_COMPONENT)
                nodal_area = node.GetSolutionStepValue(DEM.DEM_NODAL_AREA)
                total_average_target_radial_normal_stress += nodal_area * total_radial_normal_stress_value

            for element in self.dem_solution.spheres_model_part.Elements:

                total_force_X = element.GetNode(0).GetSolutionStepValue(Kratos.EXTERNAL_APPLIED_FORCE_X)
                total_force_Y = element.GetNode(0).GetSolutionStepValue(Kratos.EXTERNAL_APPLIED_FORCE_Y)
                total_force_Z = element.GetNode(0).GetSolutionStepValue(Kratos.EXTERNAL_APPLIED_FORCE_Z)
                forces_modulus = math.sqrt(total_force_X * total_force_X + total_force_Y * total_force_Y + total_force_Z * total_force_Z)
                total_forces_modulus += forces_modulus

    def DEMSolutionIsSteady(self):

        self.dem_time_steps_per_fem_time_step += 1

        if self.DEM_steps_counter == self.number_of_DEM_steps_between_steadiness_checks:

            self.DEM_steps_counter = 0

            print("\n*** DEM stationarity will now be checked...", flush=True)

            if not DEM.StationarityChecker().CheckIfVariableIsNullInModelPart(self.dem_solution.spheres_model_part, Kratos.TOTAL_FORCES_X, self.stationarity_measuring_tolerance):
                print("  F_X is larger than the desired maximum\n", flush=True)
                self.stationarity_checking_is_activated = True
                return False
            if not DEM.StationarityChecker().CheckIfVariableIsNullInModelPart(self.dem_solution.spheres_model_part, Kratos.TOTAL_FORCES_Y, self.stationarity_measuring_tolerance):
                print("  F_Y is larger than the desired maximum\n", flush=True)
                self.stationarity_checking_is_activated = True
                return False
            if not DEM.StationarityChecker().CheckIfVariableIsNullInModelPart(self.dem_solution.spheres_model_part, Kratos.TOTAL_FORCES_Z, self.stationarity_measuring_tolerance):
                print("  F_Z is larger than the desired maximum\n", flush=True)
                self.stationarity_checking_is_activated = True
                return False

            if self.stationarity_checking_is_activated:
                print("\n  ********** DEM solution is steady! ********** \n", flush=True)
                with open('poro_solution_time_vs_dem_time_vs_dems_steps.txt', 'a') as time_and_steps_file:
                    time_and_steps_file.write(str(self.poromechanics_solution.time) + " " + str(self.dem_solution.time) + " " + str(self.dem_time_steps_per_fem_time_step) + '\n')
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
