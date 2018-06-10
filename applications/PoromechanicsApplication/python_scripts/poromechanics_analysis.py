from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import time as timer
import sys
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
        print(timer.ctime())
        self.initial_time = timer.perf_counter()
        
        # Create the ModelPart
        model_part_name = parameters["problem_data"]["model_part_name"].GetString()
        self.main_model_part = Kratos.ModelPart(model_part_name)

        self.main_model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE,
                                                  parameters["problem_data"]["domain_size"].GetInt())
        self.main_model_part.ProcessInfo.SetValue(Kratos.TIME,
                                                  parameters["problem_data"]["start_time"].GetDouble())
        #TODO
        #self.main_model_part.ProcessInfo.SetValue(Kratos.DELTA_TIME, 
        #                                          parameters["problem_data"]["time_step"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosPoro.TIME_UNIT_CONVERTER, 1.0)

        super(PoromechanicsAnalysis,self).__init__(model,parameters)

        ## Import parallel modules if needed and set number of OMP threads
        parallel=Kratos.OpenMPUtils()
        parallel.SetNumThreads(parameters["problem_data"]["number_of_threads"].GetInt())
        if (self.parallel_type == "MPI"):
            import KratosMultiphysics.MetisApplication as MetisApplication
            import KratosMultiphysics.TrilinosApplication as TrilinosApplication
            print("MPI parallel configuration. OMP_NUM_THREADS =",parallel.GetNumThreads())
        else:
            print("OpenMP parallel configuration. OMP_NUM_THREADS =",parallel.GetNumThreads())
        sys.stdout.flush()

    def _CreateSolver(self):
        solver_module = __import__(self.project_parameters["solver_settings"]["solver_type"].GetString())
        solver = solver_module.CreateSolver(self.model, self.project_parameters["solver_settings"])
        return solver

    def Initialize(self):
        #TODO: seguir despues de solver
        '''
        Construct and initialize all classes and tools used in the simulation loop.
        '''

        self._SetUpRestart()

        if self.load_restart:
            self.restart_utility.LoadRestart()
        else:
            self.solver.AddVariables()
            self.solver.ImportModelPart()
            self.solver.AddDofs()

        self.model.AddModelPart(self.main_model_part)

        # Print model_part and properties
        if(self.echo_level > 1):
            print(self.main_model_part)
            for properties in self.main_model_part.Properties:
                print(properties)
        sys.stdout.flush()

        # this should let eventual derived stages modify the model after reading.
        self.ModifyInitialProperties()
        self.ModifyInitialGeometry()

        self._SetUpListOfProcesses()
        self._SetUpAnalysis()

        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

    def InitializeSolutionStep(self):

        if self.is_printing_rank:
            Kratos.Logger.PrintInfo("Poromechanics Analysis","STEP = ", self.main_model_part.ProcessInfo[Kratos.STEP])
            Kratos.Logger.PrintInfo("Poromechanics Analysis","TIME = ", self.time)

        super(PoromechanicsAnalysis,self).InitializeSolutionStep()

    def OutputSolutionStep(self):

        if self.have_output and self.output.IsOutputStep():

            for process in self.list_of_processes:
                process.ExecuteBeforeOutputStep()

            self.output.PrintOutput()

            for process in self.list_of_processes:
                process.ExecuteAfterOutputStep()

        if self.save_restart:
            self.restart_utility.SaveRestart()

    def Finalize(self):

        super(PoromechanicsAnalysis,self).Finalize()

        # Finalizing strategy
        if self.parallel_type == "OpenMP":
            self.solver.Clear()

        # Time control
        print("Analysis Completed. Elapsed Time = %.3f" % (timer.perf_counter() - self.initial_time)," seconds.")
        print(timer.ctime())

    def _SetUpListOfProcesses(self):
        '''
        Read the definition of initial and boundary conditions for the problem and initialize the processes that will manage them.
        Also initialize any additional processes present in the problem (such as those used to calculate additional results).
        '''
        from process_factory import KratosProcessFactory
        factory = KratosProcessFactory(self.model)
        # The list of processes will contain a list with each individual process already constructed
        self.list_of_processes =  factory.ConstructListOfProcesses( self.project_parameters["constraints_process_list"] )
        self.list_of_processes += factory.ConstructListOfProcesses( self.project_parameters["loads_process_list"] )

        #TODO this should be generic
        # initialize GiD  I/O
        self.output = self._SetUpGiDOutput()
        if self.output is not None:
            self.list_of_processes += [self.output,]

    def _SetUpAnalysis(self):
        '''
        Initialize the Python solver and its auxiliary tools and processes.
        This function should prepare everything so that the simulation
        can start immediately after exiting it.
        '''

        for process in self.list_of_processes:
            process.ExecuteInitialize()

        self.solver.Initialize()

        ## If the echo level is high enough, write the complete list of settings used to run the simulation
        if self.is_printing_rank and self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

        ## Stepping and time settings
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

        if self.main_model_part.ProcessInfo[Kratos.IS_RESTARTED]:
            self.time = self.main_model_part.ProcessInfo[Kratos.TIME]
        else:
            self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()


    def _SetUpGiDOutput(self):
        '''Initialize self.output as a GiD output instance.'''
        self.have_output = self.project_parameters.Has("output_configuration")
        if self.have_output:
            if self.parallel_type == "OpenMP":
                import poromechanics_cleaning_utility
                poromechanics_cleaning_utility.CleanPreviousFiles(os.getcwd()) # Clean previous post files
                from gid_output_process import GiDOutputProcess as OutputProcess
            elif self.parallel_type == "MPI":
                from gid_output_process_mpi import GiDOutputProcessMPI as OutputProcess

            output = OutputProcess(self.solver.GetComputingModelPart(),
                                   self.project_parameters["problem_data"]["problem_name"].GetString() ,
                                   self.project_parameters["output_configuration"])

            return output

    def _SetUpRestart(self):
        """Initialize self.restart_utility as a RestartUtility instance and check if we need to initialize the problem from a restart file."""
        self.load_restart = False
        self.save_restart = False

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
