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
        if (parameters["problem_data"]["parallel_type"].GetString() == "MPI"):
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

    def Initialize(self):
        if (self.reset_displacements):
            self._GetSolver().main_model_part.ProcessInfo[KratosGeo.RESET_DISPLACEMENTS] = True

        super().Initialize()
        if (self.reset_displacements):
            KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.DISPLACEMENT,self._GetSolver().GetComputingModelPart().Nodes)
            KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.ROTATION,self._GetSolver().GetComputingModelPart().Nodes)

            for node in self._GetSolver().GetComputingModelPart().Nodes:
                dNew = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 1, dNew)
                rotNew = node.GetSolutionStepValue(KratosMultiphysics.ROTATION,0)
                node.SetSolutionStepValue(KratosMultiphysics.ROTATION, 1, rotNew)

            KratosMultiphysics.VariableUtils().UpdateCurrentToInitialConfiguration(self._GetSolver().GetComputingModelPart().Nodes)

        else:
            self._GetSolver().main_model_part.ProcessInfo[KratosGeo.RESET_DISPLACEMENTS] = False

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
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "STEP: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "TIME: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME])

class GeoMechanicsAnalysis(GeoMechanicsAnalysisBase):

    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)

        self.reduction_factor    = project_parameters["solver_settings"]["reduction_factor"].GetDouble()
        self.increase_factor     = project_parameters["solver_settings"]["increase_factor"].GetDouble()
        self.delta_time          = project_parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()

        self.max_delta_time_factor = 1000.0
        if project_parameters["solver_settings"]["time_stepping"].Has("max_delta_time_factor"):
            self.max_delta_time_factor = \
                project_parameters["solver_settings"]["time_stepping"]["max_delta_time_factor"].GetDouble()

        self.max_delta_time      = self.delta_time * self.max_delta_time_factor
        self.max_iterations      = project_parameters["solver_settings"]["max_iterations"].GetInt()
        self.min_iterations      = project_parameters["solver_settings"]["min_iterations"].GetInt()
        self.number_cycles       = project_parameters["solver_settings"]["number_cycles"].GetInt()
        self.solution_type       = project_parameters["solver_settings"]["solution_type"].GetString()
        self.reset_displacements = project_parameters["solver_settings"]["reset_displacements"].GetBool()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        if(self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] > self.max_iterations):
            raise Exception("max_number_of_iterations_exceeded")

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        if self._GetSolver().settings["reset_displacements"].GetBool():
            old_total_displacements = [node.GetSolutionStepValue(KratosGeo.TOTAL_DISPLACEMENT)
                                       for node in self._GetSolver().GetComputingModelPart().Nodes]

        while self.KeepAdvancingSolutionLoop():
            if(self.delta_time > self.max_delta_time):
                self.delta_time = self.max_delta_time
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "reducing delta_time to max_delta_time: ", self.max_delta_time)
            t = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
            new_time = t + self.delta_time
            if (new_time > self.end_time):
                new_time = self.end_time
                self.delta_time = new_time - t
            self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP] += 1
            self._GetSolver().main_model_part.CloneTimeStep(new_time)
            self._GetSolver().main_model_part.ProcessInfo[KratosMultiphysics.START_TIME] = self.time
            self._GetSolver().main_model_part.ProcessInfo[KratosMultiphysics.END_TIME] = self.end_time

            converged = False
            number_cycle = 0
            while (not converged and number_cycle < self.number_cycles):

                number_cycle +=1
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "cycle: ", number_cycle)
                t = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
                corrected_time = t - self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.DELTA_TIME] + self.delta_time
                self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = corrected_time
                self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.DELTA_TIME] = self.delta_time

                self.InitializeSolutionStep()
                self._GetSolver().Predict()
                converged = self._GetSolver().SolveSolutionStep()

                if (self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] >= self.max_iterations or not converged):
                    KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Down-scaling with factor: ", self.reduction_factor)
                    self.delta_time *= self.reduction_factor

                    #converged = False
                    # Reset displacements to the initial
                    KratosMultiphysics.VariableUtils().UpdateCurrentPosition(self._GetSolver().GetComputingModelPart().Nodes, KratosMultiphysics.DISPLACEMENT,1)
                    for node in self._GetSolver().GetComputingModelPart().Nodes:
                        dold = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,1)
                        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,dold)

                    # for node in self._GetSolver().GetComputingModelPart().Nodes:
                    #     # adding TOTAL_DISPLACEMENT as dofs
                    #     KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "DISPLACEMENT_0: ", node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT))
                    #     KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "DISPLACEMENT_1: ", node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,1))

                elif (self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] < self.min_iterations):
                    KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Up-scaling with factor: ", self.increase_factor)
                    #converged = True
                    self.delta_time *= self.increase_factor
                    t = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
                    corrected_time = t + self.delta_time
                    if (corrected_time > self.end_time):
                        corrected_time = self.end_time
                        self.delta_time = corrected_time - t

            if self._GetSolver().settings["reset_displacements"].GetBool() and converged:
                for idx, node in enumerate(self._GetSolver().GetComputingModelPart().Nodes):
                    self._CalculateTotalDisplacement(node, old_total_displacements[idx])

            if (not converged):
                raise Exception('The maximum number of cycles is reached without convergence!')

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
