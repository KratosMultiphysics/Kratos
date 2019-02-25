from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Time monitoring
import time as timer

## Importing modules -----------------------------------------------------------------------------------------

# Import system python
import os

# Import kratos core and applications
import KratosMultiphysics as Kratos
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.ConvectionDiffusionApplication
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.FluidTransportApplication as KratosFluidTransport

from analysis_stage import AnalysisStage

class FluidTransportAnalysis(AnalysisStage):
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
        super(FluidTransportAnalysis,self).__init__(model,parameters)

    def _CreateSolver(self):
        solver_module = __import__(self.project_parameters["solver_settings"]["solver_type"].GetString())
        solver = solver_module.CreateSolver(self.model, self.project_parameters["solver_settings"])
        return solver

    def _GetOrderOfProcessesInitialization(self):
        return ["constraints_process_list",
                "loads_process_list",
                "auxiliar_process_list"]

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super(FluidTransportAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["constraints_process_list", "loads_process_list","auxiliar_process_list"]
            if len(list_of_processes) == 0: # Processes are given in the old format
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Using the old way to create the processes, this will be removed!")
                from process_factory import KratosProcessFactory
                factory = KratosProcessFactory(self.model)
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        list_of_processes += factory.ConstructListOfProcesses(self.project_parameters[process_name])
            else: # Processes are given in the new format
                for process_name in processes_block_names:
                    if (self.project_parameters.Has(process_name) is True):
                        raise Exception("Mixing of process initialization is not alowed!")
        elif parameter_name == "output_processes":
            if self.project_parameters.Has("output_configuration"):
                gid_output= self._SetUpGiDOutput()
                list_of_processes += [gid_output,]
        else:
            raise NameError("wrong parameter name")

        return list_of_processes

    def _SetUpGiDOutput(self):
        '''Initialize a GiD output instance.'''
        if self.parallel_type == "OpenMP":
            import fluid_transport_cleaning_utility
            fluid_transport_cleaning_utility.CleanPreviousFiles(os.getcwd()) # Clean previous post files
            from gid_output_process import GiDOutputProcess as OutputProcess
        elif self.parallel_type == "MPI":
            from gid_output_process_mpi import GiDOutputProcessMPI as OutputProcess

        output = OutputProcess(self._GetSolver().GetComputingModelPart(),
                                self.project_parameters["problem_data"]["problem_name"].GetString() ,
                                self.project_parameters["output_configuration"])

        return output

    def _GetSimulationName(self):
        return "Fluid Transport Analysis"

    def Finalize(self):
        super(FluidTransportAnalysis,self).Finalize()

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
        err_msg += '    "python fluid_transport_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python fluid_transport_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = FluidTransportAnalysis(model,parameters)
    simulation.Run()
