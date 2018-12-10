from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import time as timer
import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.PoromechanicsApplication as KratosPoro

from analysis_stage import AnalysisStage

class PoromechanicsAnalysis(AnalysisStage):
    '''Main script for poromechanics simulations.'''

    def __init__(self,model,parameters):
        # Time monitoring
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(),timer.ctime())
        self.initial_time = timer.perf_counter()

        # Set number of OMP threads
        parallel=Kratos.OpenMPUtils()
        parallel.SetNumThreads(parameters["problem_data"]["number_of_threads"].GetInt())

        ## Import parallel modules if needed
        if (parameters["problem_data"]["parallel_type"].GetString() == "MPI"):
            import KratosMultiphysics.MetisApplication as MetisApplication
            import KratosMultiphysics.TrilinosApplication as TrilinosApplication
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(),"MPI parallel configuration. OMP_NUM_THREADS =",parallel.GetNumThreads())
        else:
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(),"OpenMP parallel configuration. OMP_NUM_THREADS =",parallel.GetNumThreads())

        # Creating solver and model part and adding variables
        super(PoromechanicsAnalysis,self).__init__(model,parameters)

    def _CreateSolver(self):
        solver_module = __import__(self.project_parameters["solver_settings"]["solver_type"].GetString())
        solver = solver_module.CreateSolver(self.model, self.project_parameters["solver_settings"])
        return solver

    def _GetOrderOfProcessesInitialization(self):
        return ["constraints_process_list",
                "loads_process_list",
                "auxiliar_process_list"]

    def _GetSimulationName(self):
        return "Poromechanics Analysis"

    def Finalize(self):
        super(PoromechanicsAnalysis,self).Finalize()

        # Finalizing strategy
        if self.parallel_type == "OpenMP":
            self._GetSolver().Clear()

        # Time control
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(),"Analysis Completed. Elapsed Time = %.3f" % (timer.perf_counter() - self.initial_time)," seconds.")
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(),timer.ctime())


if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python poromechanics_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python poromechanics_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = PoromechanicsAnalysis(model,parameters)
    simulation.Run()
