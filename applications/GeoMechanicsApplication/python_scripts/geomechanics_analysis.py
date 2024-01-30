import time as timer
import os
import sys

sys.path.append(os.path.join('..','..','..'))

import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.GeoMechanicsApplication import geomechanics_solvers_wrapper

from importlib import import_module

class GeoMechanicsAnalysisBase(AnalysisStage):
    '''Main script for geomechanics simulations.'''

    def __init__(self,model,parameters):
        # Time monitoring
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(),timer.ctime())
        self.initial_time = timer.perf_counter()

        # Set number of OMP threads
        parallel=Kratos.OpenMPUtils()
        problem_data_settings = parameters["problem_data"]
        if problem_data_settings.Has("number_of_threads"):
          parallel.SetNumThreads(parameters["problem_data"]["number_of_threads"].GetInt())
        else:
          parallel.SetNumThreads(1)

        ## Import parallel modules if needed
        if parameters["problem_data"]["parallel_type"].GetString() == "MPI":
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(),"MPI parallel configuration. OMP_NUM_THREADS =",parallel.GetNumThreads())
        else:
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(),"OpenMP parallel configuration. OMP_NUM_THREADS =",parallel.GetNumThreads())

        # Creating solver and model part and adding variables
        super().__init__(model,parameters)

    def _CreateSolver(self):
        solver = geomechanics_solvers_wrapper.CreateSolver(self.model, self.project_parameters)
        return solver

    def _GetOrderOfProcessesInitialization(self):
        return ["constraints_process_list",
                "loads_process_list",
                "auxiliar_process_list"]

    def _GetSimulationName(self):
        return "GeoMechanics Analysis"

    def _CalculateTotalDisplacement(self,node, old_total_displacement):
        """
        Calculates total displacement
        :param node:
        :return:
        """
        stage_displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        total_displacement = old_total_displacement + stage_displacement
        node.SetSolutionStepValue(KratosGeo.TOTAL_DISPLACEMENT, total_displacement)

    def ResetIfHasNodalSolutionStepVariable(self, variable):
        if self._GetSolver().main_model_part.HasNodalSolutionStepVariable(variable):
            KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(variable, self._GetSolver().GetComputingModelPart().Nodes)
            for node in self._GetSolver().GetComputingModelPart().Nodes:
                new_value = node.GetSolutionStepValue(variable, 0)
                node.SetSolutionStepValue(variable, 1, new_value)

    def ModifyInitialGeometry(self):
        # Overrides the base class. Necessary to let reset_displacements function correctly i.c.w. prescribed displacements/rotations.
        # The reset needs to take place befor the Initialize of the processes, as these will set the Dirichlet condition.
        self._GetSolver().main_model_part.ProcessInfo[KratosGeo.RESET_DISPLACEMENTS] = self.reset_displacements
        if self.reset_displacements:
            self.ResetIfHasNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
            self.ResetIfHasNodalSolutionStepVariable(KratosMultiphysics.ROTATION)

            KratosMultiphysics.VariableUtils().UpdateCurrentToInitialConfiguration(self._GetSolver().GetComputingModelPart().Nodes)

    def Finalize(self):
        super().Finalize()

        # Finalizing strategy
        if self.parallel_type == "OpenMP":
            self._GetSolver().Clear()

        # Time control
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(),"Analysis Completed. Elapsed Time = %.3f" % (timer.perf_counter() - self.initial_time)," seconds.")
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(),timer.ctime())

    def KeepAdvancingSolutionLoop(self):
        return self._GetSolver().KeepAdvancingSolutionLoop(self.end_time)

    def PrintAnalysisStageProgressInformation(self):
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "STEP      : ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "DELTA_TIME: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.DELTA_TIME])
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "TIME      : ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME])

class GeoMechanicsAnalysis(GeoMechanicsAnalysisBase):

    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)

        # time step related stuff
        self.start_time          = project_parameters["problem_data"]["start_time"].GetDouble()
        self.end_time            = project_parameters["problem_data"]["end_time"].GetDouble()

        self.delta_time          = project_parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()
        self.delta_time          = min(self.delta_time, self.end_time - self.start_time)

        self.reduction_factor    = project_parameters["solver_settings"]["reduction_factor"].GetDouble()
        self.increase_factor     = project_parameters["solver_settings"]["increase_factor"].GetDouble()
        self.min_iterations      = project_parameters["solver_settings"]["min_iterations"].GetInt()
        self.max_delta_time_factor = 1000.0
        if project_parameters["solver_settings"]["time_stepping"].Has("max_delta_time_factor"):
            self.max_delta_time_factor = project_parameters["solver_settings"]["time_stepping"]["max_delta_time_factor"].GetDouble()
        self.max_delta_time      = self.delta_time * self.max_delta_time_factor
        self.number_cycles       = project_parameters["solver_settings"]["number_cycles"].GetInt()

        self.max_iterations      = project_parameters["solver_settings"]["max_iterations"].GetInt()
        self.solution_type       = project_parameters["solver_settings"]["solution_type"].GetString()
        self.reset_displacements = project_parameters["solver_settings"]["reset_displacements"].GetBool()
        self.rebuild_level       = project_parameters["solver_settings"]["rebuild_level"].GetInt()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] > self.max_iterations:
            raise Exception("max_number_of_iterations_exceeded")

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """

        # store total displacement field for reset_displacements
        if self._GetSolver().settings["reset_displacements"].GetBool():
            old_total_displacements = [node.GetSolutionStepValue(KratosGeo.TOTAL_DISPLACEMENT)
                                       for node in self._GetSolver().GetComputingModelPart().Nodes]

        self._GetSolver().solver.SetRebuildLevel(self.rebuild_level)

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
            small_time_increment = 1.E-3 * self.delta_time
            if self.end_time - new_time < small_time_increment:
                new_time = self.end_time
                self.delta_time = new_time - t
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Up-scaling to reach end_time without small increments: ", self.delta_time)

            # start the new step
            self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP] += 1
            self._GetSolver().main_model_part.CloneTimeStep(new_time)

            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "--------------------------------------", " ")

            converged = False
            number_cycle = 0
            while (not converged and number_cycle < self.number_cycles):

                number_cycle += 1
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "cycle: ", number_cycle)

                # set new_time and delta_time in the nonlinear solver
                new_time = t + self.delta_time
                self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]       = new_time
                self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.DELTA_TIME] = self.delta_time

                # do the nonlinear solver iterations
                self.InitializeSolutionStep()
                self._GetSolver().Predict()
                converged = self._GetSolver().SolveSolutionStep()
                self._GetSolver().solver.SetStiffnessMatrixIsBuilt(True)

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
                    # Reset displacements to the initial
                    KratosMultiphysics.VariableUtils().UpdateCurrentPosition(self._GetSolver().GetComputingModelPart().Nodes, KratosMultiphysics.DISPLACEMENT,1)
                    for node in self._GetSolver().GetComputingModelPart().Nodes:
                        dold = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,1)
                        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, dold)

            if not converged:
                raise Exception('The maximum number of cycles is reached without convergence!')

            if self._GetSolver().settings["reset_displacements"].GetBool():
                for idx, node in enumerate(self._GetSolver().GetComputingModelPart().Nodes):
                    self._CalculateTotalDisplacement(node, old_total_displacements[idx])

            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python geomechanics_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python geomechanics_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = GeoMechanicsAnalysis(model,parameters)
    simulation.Run()
