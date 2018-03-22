from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication
try:
    import KratosMultiphysics.ExternalSolversApplication
except ImportError:
    pass

class FluidDynamicsAnalysis(object):
    '''Main script for fluid dynamics simulations using the navier_stokes family of python solvers.'''

    def __init__(self,parameters):
        super(FluidDynamicsAnalysis,self).__init__()

        self.project_parameters = parameters

        self.echo_level = self.project_parameters["problem_data"]["echo_level"].GetInt()
        self.parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()

        # If this is an MPI run, load the distributed memory modules
        if (self.parallel_type == "MPI"):
            from KratosMultiphysics.mpi import mpi
            import KratosMultiphysics.MetisApplication
            import KratosMultiphysics.TrilinosApplication
            self.is_printing_rank = (mpi.rank == 0)
        else:
            self.is_printing_rank = True

    def SetUpModel(self):
        '''Initialize the model part for the problem and other general model data.'''

        model_part_name = self.project_parameters["problem_data"]["model_part_name"].GetString()
        self.main_model_part = Kratos.ModelPart(model_part_name)

        domain_size = self.project_parameters["problem_data"]["domain_size"].GetInt()
        self.main_model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, domain_size)

        ## Solver construction
        import python_solvers_wrapper_fluid
        self.solver = python_solvers_wrapper_fluid.CreateSolver(self.main_model_part, self.project_parameters)

        self._SetUpRestart()

        if self.load_restart:
            self.restart_utility.LoadRestart()
        else:
            self.solver.AddVariables()
            self.solver.ImportModelPart()
            self.solver.AddDofs()

        # Fill a Model instance using input
        self.model = Kratos.Model()
        self.model.AddModelPart(self.main_model_part)

    def SetUpAuxiliaryProcesses(self):
        '''
        Read the definition of initial and boundary conditions for the problem and initialize the processes that will manage them.
        Also initialize any additional processes present in the problem (such as those used to calculate additional results).
        '''
        from process_factory import KratosProcessFactory
        factory = KratosProcessFactory(self.model)
        # The list of processes will contain a list with each individual process already constructed (boundary conditions, initial conditions and gravity)
        # Note 1: gravity is constructed first. Outlet process might need its information.
        # Note 2: initial conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
        self.simulation_processes =  factory.ConstructListOfProcesses( self.project_parameters["gravity"] )
        self.simulation_processes += factory.ConstructListOfProcesses( self.project_parameters["initial_conditions_process_list"] )
        self.simulation_processes += factory.ConstructListOfProcesses( self.project_parameters["boundary_conditions_process_list"] )
        self.simulation_processes += factory.ConstructListOfProcesses( self.project_parameters["auxiliar_process_list"] )

    def SetUpAnalysis(self):
        '''
        Initialize the Python solver and its auxiliary tools and processes.
        This function should prepare everything so that the simulation
        can start immediately after exiting it.
        '''

        for process in self.simulation_processes:
            process.ExecuteInitialize()

        self.solver.Initialize()

        #TODO this should be generic
        # initialize GiD  I/O
        self._SetUpGiDOutput()

        ## Writing the full ProjectParameters file before solving
        if self.is_printing_rank and self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

        ## Stepping and time settings
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

        if self.main_model_part.ProcessInfo[Kratos.IS_RESTARTED]:
            self.time = self.main_model_part.ProcessInfo[Kratos.TIME]
            self.step = self.main_model_part.ProcessInfo[Kratos.STEP]
        else:
            self.time = 0.0
            self.step = 0

        for process in self.simulation_processes:
            process.ExecuteBeforeSolutionLoop()

        if self.have_output:
            self.output.ExecuteBeforeSolutionLoop()


    def _SetUpGiDOutput(self):
        '''Initialize self.output as a GiD output instance.'''
        self.have_output = self.project_parameters.Has("output_configuration")
        if self.have_output:
            if self.parallel_type == "OpenMP":
                from gid_output_process import GiDOutputProcess as OutputProcess
            elif self.parallel_type == "MPI":
                from gid_output_process_mpi import GiDOutputProcessMPI as OutputProcess

            self.output = OutputProcess(self.solver.GetComputingModelPart(),
                                        self.project_parameters["problem_data"]["problem_name"].GetString() ,
                                        self.project_parameters["output_configuration"])

            self.output.ExecuteInitialize()

    def _SetUpRestart(self):
        """Initialize self.restart_utility as a RestartUtility instance and check if we need to initialize the problem from a restart file."""
        if self.project_parameters.Has("restart_settings"):
            restart_settings = self.project_parameters["restart_settings"]
            self.load_restart = restart_settings["load_restart"].GetBool()
            self.save_restart = restart_settings["save_restart"].GetBool()
            restart_settings.RemoveValue("load_restart")
            restart_settings.RemoveValue("save_restart")
            restart_settings.AddValue("input_filename", self.project_parameters["problem_data"]["problem_name"])
            restart_settings.AddValue("echo_level", self.project_parameters["problem_data"]["echo_level"])

            if self.parallel_type == "OpenMP":
                from restart_utility import RestartUtility as Restart
            elif self.parallel_type == "MPI":
                from trilinos_restart_utility import TrilinosRestartUtility as Restart

            self.restart_utility = Restart(self.main_model_part,
                                           self.project_parameters["restart_settings"])
        else:
            self.load_restart = False
            self.save_restart = False

    def _TimeBufferIsInitialized(self):
        # We always have one extra old step (step 0, read from input)
        return self.step + 1 >= self.solver.GetMinimumBufferSize()

    def RunMainTemporalLoop(self):
        '''The main solution loop.'''
        while self.time <= self.end_time:

            dt = self.solver.ComputeDeltaTime()
            self.time = self.time + dt
            self.step = self.step + 1

            self.main_model_part.CloneTimeStep(self.time)
            self.main_model_part.ProcessInfo[Kratos.STEP] = self.step

            if self.is_printing_rank:
                Kratos.Logger.PrintInfo("Fluid Dynamics Analysis","STEP = ", self.step)
                Kratos.Logger.PrintInfo("Fluid Dynamics Analysis","TIME = ", self.time)

            self.InitializeSolutionStep()
            self.SolveSingleStep()
            self.FinalizeSolutionStep()

    def InitializeSolutionStep(self):

        for process in self.simulation_processes:
            process.ExecuteInitializeSolutionStep()

        if self.have_output:
            self.output.ExecuteInitializeSolutionStep()

        if self._TimeBufferIsInitialized():
            self.solver.InitializeSolutionStep()

    def SolveSingleStep(self):
        if self._TimeBufferIsInitialized():
            self.solver.SolveSolutionStep()

    def FinalizeSolutionStep(self):

        if self._TimeBufferIsInitialized():
            self.solver.FinalizeSolutionStep()

        # shouldn't this go at the end of the iteration???
        for process in self.simulation_processes:
            process.ExecuteFinalizeSolutionStep()

        if self.have_output:
            self.output.ExecuteFinalizeSolutionStep()

        if self.have_output and self.output.IsOutputStep():

            for process in self.simulation_processes:
                process.ExecuteBeforeOutputStep()

            self.output.PrintOutput()

            for process in self.simulation_processes:
                process.ExecuteAfterOutputStep()

        if self.save_restart:
            self.restart_utility.SaveRestart()


    def FinalizeAnalysis(self):
        '''Finalize the simulation and close open files.'''

        for process in self.simulation_processes:
            process.ExecuteFinalize()

        if self.have_output:
            self.output.ExecuteFinalize()

    def InitializeAnalysis(self):
        '''Wrapper function comprising the definition of the model and the initialization of the problem.'''
        self.SetUpModel()
        self.SetUpAuxiliaryProcesses()
        self.SetUpAnalysis()

    def Run(self):
        '''Wrapper function for the solution.'''
        self.InitializeAnalysis()
        self.RunMainTemporalLoop()
        self.FinalizeAnalysis()

if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python fluid_dynamics_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python fluid_dynamics_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    simulation = FluidDynamicsAnalysis(parameters)
    simulation.Run()
