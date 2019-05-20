import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import collections
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return GaussSeidel(parameters)


class GaussSeidel(object):

    def __init__(self, parameters):
        self.parameters = parameters
        self.settings = parameters["settings"]
        default_settings = cs_data_structure.Parameters("""
        {
            "name" : "",
            "solver_type" : "gauss_seidel_strong",
            "echo_level" : 0,
            "start_coupling_time" : 0.0,
            "convergence_accelerators" : [],
            "participants" : [],
            "convergence_criteria" : []
        }
        """)
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.echo_level = self.settings["echo_level"].GetInt()

        # Get the participating solvers a map with their names and objects
        self.number_of_participants = self.settings['participants'].size()
        self.participating_solvers = self._CreateSolvers(self.full_settings['solvers'])
        self.solver_settings = self._GetSolverCoSimulationDetails(self.settings["participants"])

        # With this setting the coupling can start later
        self.start_coupling_time = 0.0
        if self.settings.Has("start_coupling_time"):
            self.start_coupling_time = self.settings["start_coupling_time"].GetDouble()
        if self.start_coupling_time > 0.0:
            self.coupling_started = False
        else:
            self.coupling_started = True

    def Initialize(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.Initialize()

    def Finalize(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.Finalize()

    def Predict(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.Predict()

    def InitializeSolutionStep(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.FinalizeSolutionStep()

    def OutputSolutionStep(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.OutputSolutionStep()

    def AdvanceInTime(self, current_time):
        new_time = 0.0
        for solver_name, solver in self.participating_solvers.items():
            new_time = max(solver.AdvanceInTime(current_time), new_time)

        if self.start_coupling_time > new_time:
            self.coupling_started = False
        else:
            self.coupling_started = True

        return new_time

    def SolveSolutionStep(self):
        err_msg = 'Calling "SolveSolutionStep" of the "CoSimulationBaseCouplingSolver"!\n'
        err_msg += 'This function has to be implemented in the derived class!'
        raise NotImplementedError(err_msg)

    def Check(self):
        for solver_name, solver in self.participating_solvers.items():
            solver.Check()

    def PrintInfo(self):
        cs_tools.PrintInfo("The class", self.__class__.__name__," has the following participants")
        for solver_name, solver in self.participating_solvers.items():
            solver.PrintInfo()

    ## _SynchronizeInputData : Protected Function to obtain new (if any) input data for the
    #                          participating solvers. That will get the data necessary from
    #                          python cosim solver with name solver_name and export (if necessary)
    #                          it to the remote solver of solver_name
    #
    #
    #  So here the _SynchronizeXXXXXData functions exchange data between the python co simulation solvers.
    #  How the solvers get their data from their remote solvers is handled in ImportCouplingInterfaceData and ExportCouplingInterfaceData
    #  functions. Here Import and Export data functions of to_solver and from_solver are called.
    #
    #  @param solver_name     string: name of the solver for which data has to be synchronized
    def _SynchronizeInputData(self, solver_name):
        if self.coupling_started:
            solver = self.participating_solvers[solver_name]
            input_data_list = self.solver_settings[solver_name]["input_data_list"]
            num_input_data = input_data_list.size()

            for i in range(num_input_data):
                input_data = input_data_list[i]
                from_solver = self.participating_solvers[input_data["from_solver"].GetString()]
                data_name = input_data["data_name"].GetString()
                from_solver_data = from_solver.GetInterfaceData(data_name)
                solver_data = solver.GetInterfaceData(input_data["destination_data"].GetString())

                if(input_data.Has("mapper_settings")):
                    solver_data.origin_data = from_solver_data
                    solver_data.mapper_settings = input_data["mapper_settings"]

                solver.ImportCouplingInterfaceData(solver_data, from_solver)

    ## _SynchronizeOutputData : Protected Function to synchronize the out put data between the solver
    #                           interface and the remote solver. This assumes that the remote solver
    #                           has output the data in the format specified in the settings (of the out put data def in JSON)
    #
    def _SynchronizeOutputData(self, solver_name):
        if self.coupling_started:
            solver = self.participating_solvers[solver_name]
            output_data_list = self.solver_settings[solver_name]["output_data_list"]
            num_output_data = output_data_list.size()

            for i in range(num_output_data):
                output_data = output_data_list[i]
                to_solver = self.participating_solvers[output_data["to_solver"].GetString()]
                data_name = output_data["data_name"].GetString()
                to_solver_data = to_solver.GetInterfaceData(data_name)
                solver_data = solver.GetInterfaceData(output_data["origin_data"].GetString())

                if(output_data.Has("mapper_settings")):
                    solver_data.destination_data = to_solver_data
                    solver_data.mapper_settings = output_data["mapper_settings"]

                solver.ExportCouplingInterfaceData(solver_data, to_solver)

    def _CreateSolvers(self, SolversDataMap):
        solvers_map = collections.OrderedDict()
        import KratosMultiphysics.CoSimulationApplication.custom_solver_interfaces.co_simulation_solver_factory as cs_solver_factory

        for solver_name, settings in SolversDataMap.items():
            solver = cs_solver_factory.CreateSolverInterface(solver_name, cs_data_structure.Parameters(settings))
            solvers_map[solver_name] = solver

        return solvers_map

    def _GetSolverCoSimulationDetails(self, co_simulation_solver_settings):
        num_solvers = co_simulation_solver_settings.size()
        solver_cs_details = {}
        for i_solver in range(num_solvers):
            solver_name = co_simulation_solver_settings[i_solver]["name"].GetString()
            solver_cs_details[solver_name] = co_simulation_solver_settings[i_solver]
        return solver_cs_details

    ## _CreateFilters : Protected Function to make filter objects list and store in the datafield
    #
    #  @param conv_acc_settings dict: setting of the convergence accelerator to be make
    def _CreateFilters(self, co_simulation_solver_settings): # probably better in some utils file
        import KratosMultiphysics.CoSimulationApplication.custom_convergence_accelerators.co_simulation_convergence_accelerator_factory as cs_convergence_accelerator_factory
        num_solvers = co_simulation_solver_settings.size()
        for i_solver in range(num_solvers):
            solver_name = co_simulation_solver_settings[i_solver]["name"].GetString()
            solver = self.participating_solvers[solver_name]

            # For all the input data
            input_data_list = self.solver_settings[solver_name]["input_data_list"]
            num_input_data = input_data_list.size()
            for i in range(num_input_data):
                input_data = input_data_list[i]
                filters_list = input_data["filters"]
                destination_data = solver.GetInterfaceData(input_data["destination_data"].GetString())
                for filter in filters_list:
                    accelerator = cs_convergence_accelerator_factory.CreateConvergenceAccelerator(filter, destination_data) ## TODO: should change to filter
                    destination_data.filters.append(accelerator)

            # For all the output data
            output_data_list = self.solver_settings[solver_name]["output_data_list"]
            num_output_data = output_data_list.size()
            for i in range(num_output_data):
                output_data = output_data_list[i]
                filters_list = output_data["filters"]
                origin_data = solver.GetInterfaceData(output_data["origin_data"].GetString())
                for filter in filters_list:
                    accelerator = cs_convergence_accelerator_factory.CreateConvergenceAccelerator(filter, destination_data) ## TODO: should change to filter
                    origin_data.filters.append(accelerator)

    ## _CreateConvergenceCriteria : Private Function to make convergence criteria objects list #Comment protected
    #
    #  @param conv_acc_settings dict: setting of the convergence criteria to be make
    def _CreateConvergenceCriteria(self, conv_criteria_settings): # probably better in some utils file
        conv_criteria = []
        import KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_convergence_criteria as criteria
        num_criteria = conv_criteria_settings.size()
        for i in range(num_criteria):
            criteria_setting = conv_criteria_settings[i]
            solver_name = criteria_setting["solver"].GetString()
            solver = self.participating_solvers[solver_name]
            data_name = criteria_setting["data_name"].GetString()
            data = solver.data_map[data_name]
            criteria = criteria.Create(criteria_setting, data) # Change to use interface data
            conv_criteria.append(criteria)

        return conv_criteria