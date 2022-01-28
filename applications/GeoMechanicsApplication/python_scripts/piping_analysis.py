import time as timer
import os
import sys

sys.path.append(os.path.join('..','..','..'))

import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.GeoMechanicsApplication import geomechanics_solvers_wrapper
from geomechanics_analysis import GeoMechanicsAnalysis

from importlib import import_module

class PipingAnalysisBase(GeoMechanicsAnalysis):
    '''Main script for geomechanics simulations.'''

    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)


    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        if(self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] > self.max_iterations):
            raise Exception("max_number_of_iterations_exceeded")


    def __reduce_time_step(self):
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Down-scaling with factor: ",
                                            self.reduction_factor)
        self.delta_time *= self.reduction_factor

        # converged = False
        # Reset displacements to the initial
        KratosMultiphysics.VariableUtils().UpdateCurrentPosition(self._GetSolver().GetComputingModelPart().Nodes,
                                                                 KratosMultiphysics.DISPLACEMENT, 1)
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            dold = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 1)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, dold)

    def __increase_time_step(self):
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Up-scaling with factor: ", self.increase_factor)
        # converged = True
        self.delta_time *= self.increase_factor
        t = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        new_time = t + self.delta_time
        if (new_time > self.end_time):
            new_time = self.end_time
            self.delta_time = new_time - t

    def run_flow_calculation(self):
        converged = False
        number_cycle = 0

        while (not converged and number_cycle < self.number_cycles):
            # set solver time
            number_cycle += 1
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "cycle: ", number_cycle)
            t = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
            new_time = t - self._GetSolver().GetComputingModelPart().ProcessInfo[
                KratosMultiphysics.DELTA_TIME] + self.delta_time
            self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = new_time
            self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.DELTA_TIME] = self.delta_time

            # run solver
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            converged = self._GetSolver().SolveSolutionStep()

            # reduce time step size if no converge is reached, increase time step size if n_iterations < min iterations
            if (self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER]
                    >= self.max_iterations or not converged):
                self.__reduce_time_step()
            elif (self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER]
                  < self.min_iterations):
                self.__increase_time_step()

        if (not converged):
            raise Exception('The maximum number of cycles is reached without convergence!')

        self.FinalizeSolutionStep()
        self.OutputSolutionStep()

    def update_pipe_heights(self):
        pass

    def initialise_grow_step(self):
        pass

    def finalise_grow_step(self):
        # check if pipe should grow

        # save new pipe heights if pipe should grow

        # reset pipe heights if no piping converge is found
        self.run_flow_calculation()


    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """

        while self.KeepAdvancingSolutionLoop():
            if (self.delta_time > self.max_delta_time):
                self.delta_time = self.max_delta_time
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "reducing delta_time to max_delta_time: ", self.max_delta_time)
            t = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
            new_time = t + self.delta_time
            if (new_time > self.end_time):
                new_time = self.end_time
                self.delta_time = new_time - t
            self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP] += 1
            self._GetSolver().main_model_part.CloneTimeStep(new_time)

            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "--------------------------------------", " ")

            is_piping = True
            if is_piping:

                # todo set this as solver input
                max_piping_iterations = 500

                grow_pipe = True
                while grow_pipe:

                    self.initialise_grow_step()
                    piping_converged = False
                    piping_iter = 0

                    while not piping_converged and piping_iter < max_piping_iterations:
                        self.run_flow_calculation()

                        # update all pipe heights and check for piping convergence and check if
                        self.update_pipe_heights()

                    self.finalise_grow_step()






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
    simulation = PipingAnalysisBase(model,parameters)
    simulation.Run()
