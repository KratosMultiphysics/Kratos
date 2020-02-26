from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.FSIApplication import python_solvers_wrapper_fsi

class FSIAnalysis(AnalysisStage):
    '''Main script for FSI simulations using the FSI family of python solvers.'''

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

        # Initialize the user-provided processes
        self._AnalysisStage__CreateListOfProcesses() # Why name mangling?
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()

        self._GetSolver().Initialize()
        self.Check()

        self.ModifyAfterSolverInitialize()

        for process in self._GetListOfProcesses():
            process.ExecuteBeforeSolutionLoop()

        ## Stepping and time settings
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

        fluid_is_restarted = self.fluid_main_model_part.ProcessInfo[Kratos.IS_RESTARTED]
        structure_is_restarted = self.structure_main_model_part.ProcessInfo[Kratos.IS_RESTARTED]
        is_restarted = fluid_is_restarted and structure_is_restarted

        if is_restarted:
            fluid_time = self.fluid_main_model_part.ProcessInfo[Kratos.TIME]
            structure_time = self.structure_main_model_part.ProcessInfo[Kratos.TIME]
            if (abs(fluid_time - structure_time) > 1.0e-12):
                err_msg = 'Fluid restarting time is:' + str(fluid_time) + '\n'
                err_msg += 'Structure restarting time is:' + \
                    str(structure_time) + '\n'
                err_msg += 'Restarting time must coincide between subdomains.\n'
                raise Exception(err_msg)
            self.time = fluid_time
        else:
            self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()

        ## If the echo level is high enough, print the complete list of settings used to run the simulation
        if self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

        Kratos.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START-")

    def InitializeSolutionStep(self):

        # Since no substepping is implemented yet, structure and fluid step/time must match
        step_fluid = self.fluid_main_model_part.ProcessInfo[Kratos.STEP]
        step_structure = self.structure_main_model_part.ProcessInfo[Kratos.STEP]
        if step_fluid != step_structure:
            err_msg =  'Fluid step: '+ str(step_fluid) + '\n'
            err_msg += 'Structure step: '+ str(step_structure) + '\n'
            err_msg += 'No substepping has been implemented yet. Fluid and structure step must match.'
            raise Exception(err_msg)

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
        return python_solvers_wrapper_fsi.CreateSolver(self.model, self.project_parameters)

    def _GetSimulationName(self):
        return self.project_parameters["problem_data"]["problem_name"].GetString()

    def _GetOrderOfProcessesInitialization(self):
        return ["structure_constraints_process_list",
                "structure_contact_process_list",
                "structure_loads_process_list",
                "fluid_gravity",
                "fluid_initial_conditions_process_list",
                "fluid_boundary_conditions_process_list",
                "fluid_auxiliar_process_list"]

    def _GetOrderOfOutputProcessesInitialization(self):
        return ["gid_output"]

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
                from KratosMultiphysics.mpi.distributed_restart_utility import DistributedRestartUtility as Restart
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
                from KratosMultiphysics.mpi.distributed_restart_utility import DistributedRestartUtility as Restart

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
