import time as timer
import os
import sys

sys.path.append(os.path.join('..','..','..'))

import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.GeoMechanicsApplication import geomechanics_solvers_wrapper


class GeoMechanicsAnalysis(AnalysisStage):
    def __init__(self, model, project_parameters):
        # Time monitoring
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(),timer.ctime())
        self.initial_time = timer.perf_counter()

        problem_data_settings = project_parameters["problem_data"]
        number_of_threads = problem_data_settings["number_of_threads"].GetInt() if problem_data_settings.Has("number_of_threads") else 1
        Kratos.ParallelUtilities.SetNumThreads(number_of_threads)

        parallel_configuration = problem_data_settings["parallel_type"].GetString()
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), f"{parallel_configuration} parallel configuration | number of threads = {Kratos.ParallelUtilities.GetNumThreads()}")

        # Creating solver and model part and adding variables
        super().__init__(model, project_parameters)

        # time step related stuff
        self.start_time          = problem_data_settings["start_time"].GetDouble()
        self.end_time            = problem_data_settings["end_time"].GetDouble()

        solver_settings = project_parameters["solver_settings"]
        self.delta_time          = min(solver_settings["time_stepping"]["time_step"].GetDouble(), self.end_time - self.start_time)
        self.initial_delta_time  = self.delta_time
        self.reduction_factor    = solver_settings["reduction_factor"].GetDouble()
        self.increase_factor     = solver_settings["increase_factor"].GetDouble()
        self.min_iterations      = solver_settings["min_iterations"].GetInt()
        self.max_delta_time_factor = solver_settings["time_stepping"]["max_delta_time_factor"].GetDouble() if solver_settings["time_stepping"].Has("max_delta_time_factor") else 1000.0
        self.max_delta_time      = self.delta_time * self.max_delta_time_factor
        self.min_delta_time      = solver_settings["time_stepping"]["minimum_allowable_value"].GetDouble() if solver_settings["time_stepping"].Has("minimum_allowable_value") else None
        self.number_cycles       = solver_settings["number_cycles"].GetInt()
        self.max_iterations      = solver_settings["max_iterations"].GetInt()
        self.solution_type       = solver_settings["solution_type"].GetString()
        self.reset_displacements = solver_settings["reset_displacements"].GetBool()
        self.rebuild_level       = solver_settings["rebuild_level"].GetInt()

    def Initialize(self):
        super().Initialize()

        # Displacement and rotation variables are defined as stage displacement and rotation,
        # so they need to be reset at the start of a stage
        self.ResetIfHasNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.ResetIfHasNodalSolutionStepVariable(KratosMultiphysics.ROTATION)

        self._GetSolver().main_model_part.ProcessInfo[KratosGeo.RESET_DISPLACEMENTS] = self.reset_displacements
        if self.reset_displacements:
            self.ResetIfHasNodalSolutionStepVariable(KratosGeo.TOTAL_DISPLACEMENT)
            self.ResetIfHasNodalSolutionStepVariable(KratosGeo.TOTAL_ROTATION)

            KratosMultiphysics.VariableUtils().UpdateCurrentToInitialConfiguration(self._GetSolver().GetComputingModelPart().Nodes)

    def Finalize(self):
        super().Finalize()

        # Finalizing strategy
        if self.parallel_type == "OpenMP":
            self._GetSolver().Clear()

        # Time control
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(),"Analysis Completed. Elapsed Time = %.3f" % (timer.perf_counter() - self.initial_time)," seconds.")
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(),timer.ctime())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] > self.max_iterations:
            raise RuntimeError("max_number_of_iterations_exceeded")

    def _CheckDeltaTimeSize(self):
        min_delta_time = self._GetMinDeltaTimeValueOrDefault()
        if self.delta_time < min_delta_time:
            origin_of_value = "given" if self.min_delta_time is not None else "default"
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), f"The time step {self.delta_time} is smaller than a {origin_of_value} minimum value of {min_delta_time}")
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Please check settings in Project Parameters and Materials files.")
            raise RuntimeError('The time step is too small!')

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.START_TIME] = self.start_time
        self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.END_TIME]   = self.end_time

        self._GetSolver().solving_strategy.SetRebuildLevel(self.rebuild_level)

        while self.KeepAdvancingSolutionLoop():
            # check against max_delta_time should only be necessary here when the very first increment exceeds the maximum increment.
            if self.delta_time > self.max_delta_time:
                self.delta_time = self.max_delta_time
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Reducing delta_time to max_delta_time: ", self.max_delta_time)

            # maximize delta_time to avoid exceeding the end_time
            t               = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
            self.delta_time = min(self.delta_time, self.end_time - t)
            new_time        = t + self.delta_time

            # avoid very small remaining time steps
            if self.end_time - new_time < self._GetMinDeltaTimeValueOrDefault():
                new_time = self.end_time
                self.delta_time = new_time - t
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Up-scaling to reach end_time without small increments: ", self.delta_time)

            self._CheckDeltaTimeSize()

            # start the new step
            self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP] += 1
            self._GetSolver().main_model_part.CloneTimeStep(new_time)

            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "--------------------------------------", " ")

            converged = False
            number_cycle = 0
            while not converged and number_cycle < self.number_cycles:

                number_cycle += 1
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "cycle: ", number_cycle)

                # set new_time and delta_time in the nonlinear solver
                new_time = t + self.delta_time
                self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]             = new_time
                self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.DELTA_TIME]       = self.delta_time
                self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.NUMBER_OF_CYCLES] = number_cycle

                # do the nonlinear solver iterations
                self.InitializeSolutionStep()
                converged = self._GetSolver().SolveSolutionStep()
                self._GetSolver().solving_strategy.SetStiffnessMatrixIsBuilt(True)

                if converged:
                    # scale next step if desired
                    if new_time < self.end_time:
                        if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] < self.min_iterations:
                            # scale up next step
                            self.delta_time = min(self.increase_factor * self.delta_time, self.max_delta_time)
                            if new_time + self.delta_time <= self.end_time:
                                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Up-scaling to delta time: ", self.delta_time)
                            else:
                                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Up-scaling to reach end_time: ", self.delta_time)
                                self.delta_time = self.end_time - new_time
                        elif self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] == self.max_iterations:
                            # converged, but max_iterations reached, scale down next step
                            self.delta_time *= self.reduction_factor
                            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Down-scaling with factor: ", self.reduction_factor)
                else:
                    # scale down step and restart
                    KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Down-scaling with factor: ", self.reduction_factor)
                    self.delta_time *= self.reduction_factor
                    self._CheckDeltaTimeSize()
                    # Reset displacements to the initial
                    KratosMultiphysics.VariableUtils().UpdateCurrentPosition(self._GetSolver().GetComputingModelPart().Nodes, KratosMultiphysics.DISPLACEMENT,1)
                    for node in self._GetSolver().GetComputingModelPart().Nodes:
                        dold = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,1)
                        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, dold)

            if not converged:
                raise RuntimeError('The maximum number of cycles is reached without convergence!')

            if self._GetSolver().settings["solver_type"].GetString() == "U_Pw":
                KratosGeo.CalculateIncrementalMotionProcess(
                    self._GetSolver().GetComputingModelPart(),
                    Kratos.Parameters("""{"variable_name": "DISPLACEMENT"}""")).Execute()

                KratosGeo.CalculateTotalMotionProcess(
                    self._GetSolver().GetComputingModelPart(),
                    Kratos.Parameters("""{"variable_name": "DISPLACEMENT"}""")).Execute()

                if self._GetSolver().main_model_part.HasNodalSolutionStepVariable(KratosMultiphysics.ROTATION):
                    KratosGeo.CalculateIncrementalMotionProcess(
                        self._GetSolver().GetComputingModelPart(),
                        Kratos.Parameters("""{"variable_name": "ROTATION"}""")).Execute()

                    KratosGeo.CalculateTotalMotionProcess(
                        self._GetSolver().GetComputingModelPart(),
                        Kratos.Parameters("""{"variable_name": "ROTATION"}""")).Execute()


            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def KeepAdvancingSolutionLoop(self):
        return self._GetSolver().KeepAdvancingSolutionLoop(self.end_time)

    def ResetIfHasNodalSolutionStepVariable(self, variable):
        if self._GetSolver().main_model_part.HasNodalSolutionStepVariable(variable):
            zero_vector = Kratos.Array3([0.0, 0.0, 0.0])
            KratosGeo.NodeUtilities.AssignUpdatedVectorVariableToNodes(
                self._GetSolver().GetComputingModelPart().Nodes, variable, zero_vector, 0)
            KratosGeo.NodeUtilities.AssignUpdatedVectorVariableToNodes(
                self._GetSolver().GetComputingModelPart().Nodes, variable, zero_vector, 1)

    def PrintAnalysisStageProgressInformation(self):
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "STEP      : ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "DELTA_TIME: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.DELTA_TIME])
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "TIME      : ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME])

    def _CreateSolver(self):
        return geomechanics_solvers_wrapper.CreateSolver(self.model, self.project_parameters)

    def _GetOrderOfProcessesInitialization(self):
        return ["constraints_process_list",
                "loads_process_list",
                "auxiliary_process_list"]

    def _GetSimulationName(self):
        return "GeoMechanics Analysis"

    def _GetMinDeltaTimeValueOrDefault(self):
        delta_time_as_fraction_of_time_span = 0.0001
        return self.min_delta_time if self.min_delta_time is not None else min(self.initial_delta_time, delta_time_as_fraction_of_time_span * (self.end_time - self.start_time))

if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python geomechanics_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python geomechanics_analysis.py <my-parameter-file>.json"\n'
        raise TypeError(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = GeoMechanicsAnalysis(model,parameters)
    simulation.Run()
