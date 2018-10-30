from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
import KratosMultiphysics.MeshMovingApplication as MeshMovingApplication
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
try:
    import KratosMultiphysics.ExternalSolversApplication
except ImportError:
    pass

from analysis_stage import AnalysisStage

class FSIAnalysis(AnalysisStage):
    '''Main script for FSI simulations using the FSI family of python solvers.'''

    def __init__(self, model, project_parameters):
        '''The constructor of the FSI analysis object

        Note that the base class analysis stage constructor
        is not called intentionally, since it is not compatible
        with the current FSI application Json structure.
        '''
        if (type(model) != Kratos.Model):
            raise Exception("Input is expected to be provided as a Kratos Model object")

        if (type(project_parameters) != Kratos.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.model = model
        self.project_parameters = project_parameters

        self.echo_level = self.project_parameters["solver_settings"]["echo_level"].GetInt()
        self.parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()

        # If this is an MPI run, load the distributed memory modules
        if (self.parallel_type == "MPI"):
            from KratosMultiphysics.mpi import mpi
            import KratosMultiphysics.MetisApplication
            import KratosMultiphysics.TrilinosApplication
            import KratosMultiphysics.MappingApplication #TODO: Import always once we use the serial version of the mapper
            self.is_printing_rank = (mpi.rank == 0)
        else:
            self.is_printing_rank = True

        # Deprecation warnings
        # This makes possible the FSI solver derivation from the core base python_solver.py
        # if not self.project_parameters.Has("echo_level"):
        #     Kratos.Logger.PrintInfo("FSIAnalysis", "Using the old way to pass the echo_level, this will be removed!")
        #     self.project_parameters.AddEmptyValue("echo_level")
        #     self.project_parameters["echo_level"].SetInt(self.echo_level)

        if not self.project_parameters["solver_settings"].Has("parallel_type"):
            self.project_parameters["solver_settings"].AddEmptyValue("parallel_type")
            self.project_parameters["solver_settings"]["parallel_type"].SetString(self.parallel_type)

        # Add solver variables (note that the solver is created in the first _GetSolver() call)
        self._GetSolver().AddVariables()

    def Initialize(self):
        '''
        Construct and initialize all classes and tools used in the simulation loop.
        '''
        self._SetUpRestart()

        if self.load_restart:
            fluid_restart_utility = self._GetFluidRestartUtility()
            structure_restart_utility = self._GetStructureRestartUtility()
            fluid_restart_utility.LoadRestart()
            structure_restart_utility.LoadRestart()
        else:
            self._GetSolver().ImportModelPart()
            self._GetSolver().PrepareModelPart()
            self._GetSolver().AddDofs()
            self.fluid_main_model_part = self.model.GetModelPart(
                self.project_parameters["solver_settings"]["fluid_solver_settings"]["model_part_name"].GetString())
            self.structure_main_model_part = self.model.GetModelPart(
                self.project_parameters["solver_settings"]["structure_solver_settings"]["model_part_name"].GetString())

        # This should let eventual derived stages modify the model after reading.
        self.ModifyInitialProperties()
        self.ModifyInitialGeometry()

        self._GetListOfProcesses()
        self._SetUpAnalysis()

        for process in self._GetListOfProcesses():
            process.ExecuteBeforeSolutionLoop()

    def InitializeSolutionStep(self):

        # Since no substepping is implemented yet, structure and fluid step/time must match
        step_fluid = self.fluid_main_model_part.ProcessInfo[KratosMultiphysics.STEP]
        step_structure = self.structure_main_model_part.ProcessInfo[KratosMultiphysics.STEP]
        if step_fluid != step_structure:
            err_msg =  'Fluid step: '+ str(step_fluid) + '\n'
            err_msg += 'Structure step: '+ str(step_structure) + '\n'
            err_msg += 'No substepping has been implemented yet. Fluid and structure step must match.'
            raise Exception(err_msg)

        if self.is_printing_rank:
            Kratos.Logger.PrintInfo(self._GetSimulationName(),"STEP = ", step_fluid)
            Kratos.Logger.PrintInfo(self._GetSimulationName(),"TIME = ", self.time)

        self.ApplyBoundaryConditions() #here the processes are called
        self.ChangeMaterialProperties() #this is normally empty
        self._GetSolver().InitializeSolutionStep()

    def OutputSolutionStep(self):
        super(FSIAnalysis, self).OutputSolutionStep()

        if self.save_restart:
            fluid_restart_utility = self._GetFluidRestartUtility()
            structure_restart_utility = self._GetStructureRestartUtility()
            fluid_restart_utility.SaveRestart()
            structure_restart_utility.SaveRestart()

    def _CreateSolver(self):
        import python_solvers_wrapper_fsi
        return python_solvers_wrapper_fsi.CreateSolver(self.model, self.project_parameters)

    def _GetSimulationName(self):
        return self.project_parameters["problem_data"]["problem_name"].GetString()

    def _GetOrderOfProcessesInitialization(self):
        return ["structure_constraints_process_list",
                "structure_loads_process_list",
                "fluid_gravity",
                "fluid_initial_conditions_process_list",
                "fluid_boundary_conditions_process_list",
                "fluid_auxiliar_process_list"]

    def _GetOrderOfOutputProcessesInitialization(self):
        return ["gid_output"]

    def _SetUpAnalysis(self):
        '''
        Initialize the Python solver and its auxiliary tools and processes.
        This function should prepare everything so that the simulation
        can start immediately after exiting it.
        '''

        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()

        self._GetSolver().Initialize()

        self.ModifyAfterSolverInitialize()

        ## If the echo level is high enough, print the complete list of settings used to run the simulation
        if self.is_printing_rank and self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

        ## Stepping and time settings
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

        fluid_is_restarted = self.fluid_main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]
        structure_is_restarted = self.structure_main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]
        is_restarted = fluid_is_restarted and structure_is_restarted

        if is_restarted:
            fluid_time = self.fluid_main_model_part.ProcessInfo[KratosMultiphysics.TIME]
            structure_time = self.structure_main_model_part.ProcessInfo[KratosMultiphysics.TIME]
            if (fluid_time != structure_time):
                err_msg =  'Fluid restarting time is:' + str(fluid_time) + '\n'
                err_msg += 'Structure restarting time is:' + str(structure_time) + '\n'
                err_msg += 'Restarting time must coincide between subdomains.\n'
                raise Exception(err_msg)
            self.time = fluid_time
        else:
            self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()

    def _SetUpRestart(self):
        """Initialize self.restart_utility as a RestartUtility instance and check if we need to initialize the problem from a restart file."""
        has_restart = self.project_parameters.Has("restart_settings")

        if has_restart:
            raise Exception("FSI restart not implemented yet.")
        else:
            self.load_restart = False
            self.save_restart = False

    def _GetFluidRestartUtility(self):

        if self.__fluid_restart_utility is not None:
            return self.__fluid_restart_utility
        else:
            if self.parallel_type == "OpenMP":
                from restart_utility import RestartUtility as Restart
            elif self.parallel_type == "MPI":
                from trilinos_restart_utility import TrilinosRestartUtility as Restart

            model_part_name = self.project_parameters["fluid_solver_settings"]["solver_settings"]["model_part_name"].GetString()
            if self.model.HasModelPart(model_part_name):
                model_part = self.model.GetModelPart(model_part_name)
            else:
                model_part = self.model.CreateModelPart(model_part_name)

            self.__fluid_restart_utility = Restart(
                model_part,
                self.project_parameters["fluid_solver_settings"]["restart_settings"])

    def _GetStructureRestartUtility(self):

        if self.__structure_restart_utility is not None:
            return self.__structure_restart_utility
        else:
            if self.parallel_type == "OpenMP":
                from restart_utility import RestartUtility as Restart
            elif self.parallel_type == "MPI":
                from trilinos_restart_utility import TrilinosRestartUtility as Restart

            model_part_name = self.project_parameters["structure_solver_settings"]["solver_settings"]["model_part_name"].GetString()
            if self.model.HasModelPart(model_part_name):
                model_part = self.model.GetModelPart(model_part_name)
            else:
                model_part = self.model.CreateModelPart(model_part_name)

            self.__structure_restart_utility = Restart(
                model_part,
                self.project_parameters["structure_solver_settings"]["restart_settings"])


if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python fsi_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python fsi_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = FSIAnalysis(model, parameters)
    simulation.Run()
