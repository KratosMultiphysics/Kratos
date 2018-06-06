from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KFluid
import KratosMultiphysics.HDF5Application as KHdf5
import KratosMultiphysics.AdjointFluidApplication as KAdjoint
try:
    import KratosMultiphysics.ExternalSolversApplication
except ImportError:
    pass

from analysis_stage import AnalysisStage
from fluid_dynamics_analysis import FluidDynamicsAnalysis

class AdjointFluidAnalysis(AnalysisStage):
    '''Main script for adjoint sensitivity optimization in fluid dynamics simulations.'''

    def __init__(self,model,parameters):
        super(AdjointFluidAnalysis,self).__init__(model,parameters)

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

        ## Create model part and solver (but don't initialize them yet)
        model_part_name = self.project_parameters["problem_data"]["model_part_name"].GetString()
        self.main_model_part = Kratos.ModelPart(model_part_name)

        #import python_solvers_wrapper_fluid
        #self.solver = python_solvers_wrapper_fluid.CreateSolver(self.main_model_part, self.project_parameters)
        #TODO: add adjoint solver to wrapper as part of the migration procedure
        solver_module = __import__(self.project_parameters["solver_settings"]["solver_type"].GetString())
        self.solver = solver_module.CreateSolver(self.main_model_part, self.project_parameters["solver_settings"])

    def Initialize(self):
        '''
        Construct and initialize all classes and tools used in the simulation loop.
        '''
        domain_size = self.project_parameters["problem_data"]["domain_size"].GetInt()
        self.main_model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, domain_size)

        self._SetUpRestart()

        if self.load_restart:
            self.restart_utility.LoadRestart()
        else:
            self.solver.AddVariables()
            self.solver.ImportModelPart()
            self.solver.AddDofs()

        self.model.AddModelPart(self.main_model_part)

        # this should let eventual derived stages modify the model after reading.
        self.ModifyInitialProperties()
        self.ModifyInitialGeometry()

        self._SetUpListOfProcesses()
        self._SetUpAnalysis()

        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

    def InitializeSolutionStep(self):

        if self.is_printing_rank:
            Kratos.Logger.PrintInfo("Adjoint Fluid Analysis","STEP = ", self.main_model_part.ProcessInfo[Kratos.STEP])
            Kratos.Logger.PrintInfo("Adjoint Fluid Analysis","TIME = ", self.time)

        super(AdjointFluidAnalysis,self).InitializeSolutionStep()

    def OutputSolutionStep(self):

        if self.have_output and self.output.IsOutputStep():

            for process in self.list_of_processes:
                process.ExecuteBeforeOutputStep()

            self.output.PrintOutput()

            for process in self.list_of_processes:
                process.ExecuteAfterOutputStep()

        if self.save_restart:
            self.restart_utility.SaveRestart()


    def RunSolutionLoop(self):
        """Note that the adjoint problem is solved in reverse time
        """

        #self.output.PrintOutput()

        for step in range(self.number_of_steps):
            self.time = self.solver.AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self.solver.Predict()
            self.solver.SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def _SetUpListOfProcesses(self):
        '''
        Read the definition of initial and boundary conditions for the problem and initialize the processes that will manage them.
        Also initialize any additional processes present in the problem (such as those used to calculate additional results).
        '''
        from process_factory import KratosProcessFactory
        factory = KratosProcessFactory(self.model)
        self.list_of_processes =  factory.ConstructListOfProcesses( self.project_parameters["initial_conditions_process_list"] )
        self.list_of_processes += factory.ConstructListOfProcesses( self.project_parameters["boundary_conditions_process_list"] )
        self.list_of_processes += factory.ConstructListOfProcesses( self.project_parameters["gravity"] )
        if self.project_parameters.Has("list_other_processes"):
            self.list_of_processes += factory.ConstructListOfProcesses( self.project_parameters["list_other_processes"] )

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
        self.solver.SolverInitialize() # this initializes the strategy

        ## If the echo level is high enough, print the complete list of settings used to run the simualtion
        if self.is_printing_rank and self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

        ## Stepping and time settings
        self.number_of_steps = self.project_parameters["problem_data"]["nsteps"].GetInt()

        if self.main_model_part.ProcessInfo[Kratos.IS_RESTARTED]:
            self.time = self.main_model_part.ProcessInfo[Kratos.TIME]
        else:
            self.time = self.project_parameters["problem_data"]["start_step"].GetDouble()
            self.main_model_part.ProcessInfo.SetValue(Kratos.TIME,self.time)


    def _SetUpGiDOutput(self):
        '''Initialize self.output as a GiD output instance.'''
        self.have_output = self.project_parameters.Has("output_configuration")
        if self.have_output:
            if self.parallel_type == "OpenMP":
                from gid_output_process import GiDOutputProcess as OutputProcess
            elif self.parallel_type == "MPI":
                from gid_output_process_mpi import GiDOutputProcessMPI as OutputProcess

            output = OutputProcess(self.solver.GetComputingModelPart(),
                                   self.project_parameters["problem_data"]["problem_name"].GetString() ,
                                   self.project_parameters["output_configuration"])

            return output

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

if __name__ == '__main__':
    from sys import argv

    primal_parameter_file_name = None
    adjoint_parameter_file_name = None

    parameters = Kratos.Parameters(r'''{}''')

    if len(argv) == 2:
        adjoint_parameter_file_name = argv[1]
    elif len(argv) == 3:
        primal_parameter_file_name = argv[1]
        adjoint_parameter_file_name = argv[2]
    else:
        err_msg =  'Unexpected amount of input arguments!\n'
        err_msg += 'To run the primal fluid problem followed by the adjoint solution, provide both parameter files:\n'
        err_msg += '    "python adjoint_fluid_analysis.py <primal-parameter-file>.json <adjoint-parameter-file>.json"\n'
        err_msg += 'To run only the adjoint problem, provide only the adjoint parameter file:\n'
        err_msg += '    "python adjoint_fluid_analysis.py <adjoint-parameter-file>.json"\n'
        raise Exception(err_msg)

    if primal_parameter_file_name is not None:
        with open(primal_parameter_file_name,'r') as primal_parameter_file:
            parameters.AddValue("primal_settings",Kratos.Parameters(primal_parameter_file.read()))
    else:
        parameters.AddEmptyValue("primal_settings")

    with open(adjoint_parameter_file_name,'r') as adjoint_parameter_file:
        parameters.AddValue("adjoint_settings", Kratos.Parameters(adjoint_parameter_file.read()))

    model = Kratos.Model()

    if primal_parameter_file_name is not None:
        primal_simulation = FluidDynamicsAnalysis(model,parameters["primal_settings"])
        primal_simulation.Run()

    adjoint_model = Kratos.Model()
    adjoint_simulation = AdjointFluidAnalysis(adjoint_model,parameters["adjoint_settings"])
    adjoint_simulation.Run()
