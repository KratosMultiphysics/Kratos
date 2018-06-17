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

        # Creating solver and model part and adding variables
        super(PoromechanicsAnalysis,self).__init__(model,parameters)

        ## Import parallel modules if needed
        if (self.parallel_type == "MPI"):
            import KratosMultiphysics.MetisApplication as MetisApplication
            import KratosMultiphysics.TrilinosApplication as TrilinosApplication
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(),"MPI parallel configuration. OMP_NUM_THREADS =",parallel.GetNumThreads())
        else:
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(),"OpenMP parallel configuration. OMP_NUM_THREADS =",parallel.GetNumThreads())

    def _CreateSolver(self):
        solver_module = __import__(self.project_parameters["solver_settings"]["solver_type"].GetString())
        solver = solver_module.CreateSolver(self.model, self.project_parameters["solver_settings"])
        return solver

    def _GetOrderOfProcessesInitialization(self):
        return ["constraints_process_list",
                "loads_process_list"]

    def _CreateProcesses(self, parameter_name, initialization_order):
        ### TODO: seguir...
        """Create a list of Processes
        This method is TEMPORARY to not break existing code
        It will be removed in the future
        """
        list_of_processes = super(PoromechanicsAnalysis, self)._CreateProcesses(parameter_name, initialization_order)

        if parameter_name == "processes":
            processes_block_names = ["constraints_process_list", "loads_process_list"]
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




    # def Initialize(self):
    #     '''
    #     Construct and initialize all classes and tools used in the simulation loop.
    #     '''

    #     self._SetUpRestart()

    #     if self.load_restart:
    #         self.restart_utility.LoadRestart()
    #     else:
    #         self.solver.AddVariables()
    #         self.solver.ImportModelPart()
    #         self.solver.AddDofs()

    #     self.model.AddModelPart(self.main_model_part)

    #     # Print model_part and properties
    #     if(self.echo_level > 1):
    #         print(self.main_model_part)
    #         for properties in self.main_model_part.Properties:
    #             print(properties)

    #     # this should let eventual derived stages modify the model after reading.
    #     self.ModifyInitialProperties()
    #     self.ModifyInitialGeometry()

    #     self._SetUpListOfProcesses()
    #     self._SetUpAnalysis()

    #     for process in self.list_of_processes:
    #         process.ExecuteBeforeSolutionLoop()

    # def InitializeSolutionStep(self):

    #     if self.is_printing_rank:
    #         Kratos.Logger.PrintInfo("Poromechanics Analysis","STEP = ", self.main_model_part.ProcessInfo[Kratos.STEP])
    #         Kratos.Logger.PrintInfo("Poromechanics Analysis","TIME = ", self.time)

    #     super(PoromechanicsAnalysis,self).InitializeSolutionStep()

    # def OutputSolutionStep(self):

    #     if self.have_output and self.output.IsOutputStep():

    #         for process in self.list_of_processes:
    #             process.ExecuteBeforeOutputStep()

    #         self.output.PrintOutput()

    #         for process in self.list_of_processes:
    #             process.ExecuteAfterOutputStep()

    #     if self.save_restart:
    #         self.restart_utility.SaveRestart()

    def Finalize(self):

        super(PoromechanicsAnalysis,self).Finalize()

        # Finalizing strategy
        if self.parallel_type == "OpenMP":
            self.solver.Clear()

        # Time control
        print("Analysis Completed. Elapsed Time = %.3f" % (timer.perf_counter() - self.initial_time)," seconds.")
        print(timer.ctime())

    # def _SetUpListOfProcesses(self):
    #     '''
    #     Read the definition of initial and boundary conditions for the problem and initialize the processes that will manage them.
    #     Also initialize any additional processes present in the problem (such as those used to calculate additional results).
    #     '''
    #     from process_factory import KratosProcessFactory
    #     factory = KratosProcessFactory(self.model)
    #     # The list of processes will contain a list with each individual process already constructed
    #     self.list_of_processes =  factory.ConstructListOfProcesses( self.project_parameters["constraints_process_list"] )
    #     self.list_of_processes += factory.ConstructListOfProcesses( self.project_parameters["loads_process_list"] )

    #     #TODO this should be generic
    #     # initialize GiD  I/O
    #     self.output = self._SetUpGiDOutput()
    #     if self.output is not None:
    #         self.list_of_processes += [self.output,]

    # def _SetUpAnalysis(self):
    #     '''
    #     Initialize the Python solver and its auxiliary tools and processes.
    #     This function should prepare everything so that the simulation
    #     can start immediately after exiting it.
    #     '''

    #     for process in self.list_of_processes:
    #         process.ExecuteInitialize()

    #     self.solver.Initialize()

    #     ## If the echo level is high enough, write the complete list of settings used to run the simulation
    #     if self.is_printing_rank and self.echo_level > 1:
    #         with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
    #             parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

    #     ## Stepping and time settings
    #     self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

    #     if self.main_model_part.ProcessInfo[Kratos.IS_RESTARTED]:
    #         self.time = self.main_model_part.ProcessInfo[Kratos.TIME]
    #     else:
    #         self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()


    def _SetUpGiDOutput(self):
        '''Initialize a GiD output instance.'''
        if self.parallel_type == "OpenMP":
            import poromechanics_cleaning_utility
            poromechanics_cleaning_utility.CleanPreviousFiles(os.getcwd()) # Clean previous post files
            from gid_output_process import GiDOutputProcess as OutputProcess
        elif self.parallel_type == "MPI":
            from gid_output_process_mpi import GiDOutputProcessMPI as OutputProcess

        output = OutputProcess(self._GetSolver().GetComputingModelPart(),
                                self.project_parameters["problem_data"]["problem_name"].GetString() ,
                                self.project_parameters["output_configuration"])

        return output

    # def _SetUpRestart(self):
    #     """Initialize self.restart_utility as a RestartUtility instance and check if we need to initialize the problem from a restart file."""
    #     self.load_restart = False
    #     self.save_restart = False

    def _GetSimulationName(self):
        return "Poromechanics Analysis"

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
