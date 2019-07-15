from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from . import co_simulation_solver_wrapper

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

class CoSimulationCoupledSolver(co_simulation_solver_wrapper.CoSimulationSolverWrapper):
    """Baseclass for the coupled solvers used for CoSimulation
    Performs basic operations that are common among coupled solvers:
    - holds Predictors
    - holds DataTransferOperators
    - holds CouplingOperations
    - initialization of IOs of solvers
    - Snychronization of Input and Output
    - Handles the coupling sequence
    """
    def __init__(self, settings, solver_name):
        super(CoSimulationCoupledSolver, self).__init__(settings, solver_name)

        self.solver_wrappers = self.__CreateSolverWrappers()
        self._AllocateHistoricalVariablesFromCouplingData()

        self.coupling_sequence = self.__GetSolverCoSimulationDetails()

        ### Creating the predictors
        self.predictors_list = cs_tools.CreatePredictors(
            self.settings["predictors"],
            self.solver_wrappers,
            self.echo_level)

        ### Creating the coupling operations
        self.coupling_operations_dict = cs_tools.CreateCouplingOperations(
            self.settings["coupling_operations"],
            self.solver_wrappers,
            self.echo_level)

        ### Creating the data transfer operators
        self.data_transfer_operators_dict = cs_tools.CreateDataTransferOperators(
            self.settings["data_transfer_operators"],
            self.echo_level)

    def Initialize(self):
        for solver in self.solver_wrappers.values():
            solver.Initialize()

        for solver in self.solver_wrappers.values():
            solver.InitializeIO(self.solver_wrappers, self.echo_level)
            # we use the Echo_level of the coupling solver, since IO is needed by the coupling
            # and not by the (physics-) solver

        super(CoSimulationCoupledSolver, self).Initialize()

        for predictor in self.predictors_list:
            predictor.Initialize()

        for coupling_operation in self.coupling_operations_dict.values():
            coupling_operation.Initialize()

    def Finalize(self):
        for solver in self.solver_wrappers.values():
            solver.Finalize()

        for predictor in self.predictors_list:
            predictor.Finalize()

        for coupling_operation in self.coupling_operations_dict.values():
            coupling_operation.Finalize()

    def AdvanceInTime(self, current_time):
        self.time = 0.0
        for solver in self.solver_wrappers.values():
            self.time = max(self.time, solver.AdvanceInTime(current_time))

        return self.time

    def Predict(self):
        for predictor in self.predictors_list:
            predictor.Predict()

        for solver in self.solver_wrappers.values():
            solver.Predict()

    def InitializeSolutionStep(self):
        for solver in self.solver_wrappers.values():
            solver.InitializeSolutionStep()

        for predictor in self.predictors_list:
            predictor.InitializeSolutionStep()

        for coupling_operation in self.coupling_operations_dict.values():
            coupling_operation.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        for solver in self.solver_wrappers.values():
            solver.FinalizeSolutionStep()

        for predictor in self.predictors_list:
            predictor.FinalizeSolutionStep()

        for coupling_operation in self.coupling_operations_dict.values():
            coupling_operation.FinalizeSolutionStep()

    def OutputSolutionStep(self):
        for solver in self.solver_wrappers.values():
            solver.OutputSolutionStep()

    def SolveSolutionStep(self):
        err_msg  = 'Calling "SolveSolutionStep" of the "CoSimulationCoupledSolver"!\n'
        err_msg += 'This function has to be implemented in the derived class!'
        raise Exception(err_msg)

    def _SynchronizeInputData(self, solver_name):
        to_solver = self.solver_wrappers[solver_name]
        input_data_list = self.coupling_sequence[solver_name]["input_data_list"]
        if self.echo_level > 2:
            cs_tools.cs_print_info(self._ClassName(), 'Start Synchronizing Input for "{}"'.format(colors.blue(solver_name)))

        for i in range(input_data_list.size()):
            i_input_data = input_data_list[i]

            data_name = i_input_data["data"].GetString()
            from_solver_name = i_input_data["from_solver"].GetString()
            from_solver_data_name = i_input_data["from_solver_data"].GetString()

            if self.echo_level > 2:
                cs_tools.cs_print_info("  Data", '"{}" | from solver: "{}": "{}"'.format(colors.magenta(data_name), colors.blue(from_solver_name), colors.magenta(from_solver_data_name)))

            # check if data-exchange is specified for current time
            if not KM.IntervalUtility(i_input_data).IsInInterval(self.time):
                if self.echo_level > 2:
                    cs_tools.cs_print_info("  Skipped", 'not in interval')
                continue

            # from solver
            from_solver = self.solver_wrappers[from_solver_name]
            from_solver_data = from_solver.GetInterfaceData(from_solver_data_name)

            # to solver
            to_solver_data = to_solver.GetInterfaceData(data_name)

            # Importing data from external solvers
            # to_solver.ImportCouplingInterfaceData(to_solver_data)

            # perform the data transfer
            self.__ExecuteCouplingOperations(i_input_data["before_data_transfer_operations"])

            data_transfer_operator_name = i_input_data["data_transfer_operator"].GetString()
            # TODO check the order of solvers!
            self.__GetDataTransferOperator(data_transfer_operator_name).TransferData(from_solver_data, to_solver_data, i_input_data["data_transfer_operator_options"])

            self.__ExecuteCouplingOperations(i_input_data["after_data_transfer_operations"])

            self.__ApplyScaling(to_solver_data, i_input_data)

        if self.echo_level > 2:
            cs_tools.cs_print_info(self._ClassName(), 'End Synchronizing Input for "{}"'.format(colors.blue(solver_name)))

    def _SynchronizeOutputData(self, solver_name):
        from_solver = self.solver_wrappers[solver_name]
        output_data_list = self.coupling_sequence[solver_name]["output_data_list"]
        if self.echo_level > 2:
            cs_tools.cs_print_info(self._ClassName(), 'Start Synchronizing Output for "{}"'.format(colors.blue(solver_name)))

        for i in range(output_data_list.size()):
            i_output_data = output_data_list[i]

            data_name = i_output_data["data"].GetString()
            to_solver_name = i_output_data["to_solver"].GetString()
            to_solver_data_name = i_output_data["to_solver_data"].GetString()

            if self.echo_level > 2:
                cs_tools.cs_print_info("  Data", '"{}" | to solver: "{}": "{}"'.format(colors.magenta(data_name), colors.blue(to_solver_name), colors.magenta(to_solver_data_name)))

            # check if data-exchange is specified for current time
            if not KM.IntervalUtility(i_output_data).IsInInterval(self.time):
                if self.echo_level > 2:
                    cs_tools.cs_print_info("  Skipped", 'not in interval')
                continue

            # from solver
            from_solver_data = from_solver.GetInterfaceData(data_name)

            # to solver
            to_solver = self.solver_wrappers[to_solver_name]
            to_solver_data = to_solver.GetInterfaceData(to_solver_data_name)

            # perform the data transfer
            self.__ExecuteCouplingOperations(i_output_data["before_data_transfer_operations"])

            data_transfer_operator_name = i_output_data["data_transfer_operator"].GetString()
            # TODO check the order of solvers!
            self.__GetDataTransferOperator(data_transfer_operator_name).TransferData(from_solver_data, to_solver_data, i_output_data["data_transfer_operator_options"])

            self.__ExecuteCouplingOperations(i_output_data["after_data_transfer_operations"])

            self.__ApplyScaling(to_solver_data, i_output_data)

            # Exporting data to external solvers
            # from_solver.ExportCouplingInterfaceData(from_solver_data)

        if self.echo_level > 2:
            cs_tools.cs_print_info(self._ClassName(), 'End Synchronizing Output for "{}"'.format(colors.blue(solver_name)))


    def __GetDataTransferOperator(self, data_transfer_operator_name):
        try:
            return self.data_transfer_operators_dict[data_transfer_operator_name]
        except KeyError:
            raise NameError('The data-transfer-operator "{}" does not exist!'.format(data_transfer_operator_name))


    def __ExecuteCouplingOperations(self, settings):
        for i in range(settings.size()):
            coupling_operation_name = settings[i].GetString()
            self.coupling_operations_dict[coupling_operation_name].Execute()

    def PrintInfo(self):
        super(CoSimulationCoupledSolver, self).PrintInfo()

        cs_tools.cs_print_info(self._ClassName(), "Has the following components:")
        for solver in self.solver_wrappers.values():
            solver.PrintInfo()

        for predictor in self.predictors_list:
            predictor.PrintInfo()

        for coupling_operation in self.coupling_operations_dict.values():
            coupling_operation.PrintInfo()

    def Check(self):
        super(CoSimulationCoupledSolver, self).Check()
        for solver in self.solver_wrappers.values():
            solver.Check()

        for predictor in self.predictors_list:
            predictor.Check()

        for coupling_operation in self.coupling_operations_dict.values():
            coupling_operation.Check()

    def __CreateSolverWrappers(self):
        ### ATTENTION, big flaw, also the participants can be coupled solvers !!!
        import KratosMultiphysics.CoSimulationApplication.factories.solver_wrapper_factory as solvers_wrapper_factory
        from collections import OrderedDict
        # first create all solvers
        solvers = {}
        for solver_name, solver_settings in self.settings["solvers"].items():
            solvers[solver_name] = solvers_wrapper_factory.CreateSolverWrapper(solver_settings, solver_name)

        # then order them according to the coupling-loop
        # NOTE solvers that are not used in the coupling-sequence will not participate
        solvers_map = OrderedDict()
        for i_solver_settings in range(self.settings["coupling_sequence"].size()):
            solver_settings = self.settings["coupling_sequence"][i_solver_settings]
            solver_name = solver_settings["name"].GetString()
            solvers_map[solver_name] = solvers[solver_name]

        return solvers_map

    def __GetSolverCoSimulationDetails(self):
        def ValidateAndAssignDefaultsDataList(data_list, defaults):
            for i_data_list in range(data_list.size()):
                cur_data = data_list[i_data_list]

                # doing some tricks since the type of "scaling_factor" can be double or string and hence would fail in the validation
                scaling_function_string = None
                if cur_data.Has("scaling_factor") and cur_data["scaling_factor"].IsString():
                    scaling_function_string = cur_data["scaling_factor"].GetString()
                    cur_data.RemoveValue("scaling_factor")

                cur_data.ValidateAndAssignDefaults(defaults)

                if scaling_function_string is not None:
                    cur_data["scaling_factor"].SetString(scaling_function_string)

        solver_cosim_details = {}
        for i_solver_settings in range(self.settings["coupling_sequence"].size()):
            solver_settings = self.settings["coupling_sequence"][i_solver_settings]
            solver_name = solver_settings["name"].GetString()
            solver_cosim_details[solver_name] = solver_settings

            ValidateAndAssignDefaultsDataList(solver_settings["input_data_list"], GetInputDataDefaults())
            ValidateAndAssignDefaultsDataList(solver_settings["output_data_list"], GetOutputDataDefaults())

        return solver_cosim_details

    def __ApplyScaling(self, interface_data, data_configuration):
        # perform scaling of data if specified
        if data_configuration["scaling_factor"].IsString():
            from KratosMultiphysics.CoSimulationApplication.function_callback_utility import GenericCallFunction
            scaling_function_string = data_configuration["scaling_factor"].GetString()
            scope_vars = {'t' : self.time} # make time useable in function
            scaling_factor = GenericCallFunction(scaling_function_string, scope_vars) # evaluating function string
        else:
            scaling_factor = data_configuration["scaling_factor"].GetDouble()

        if abs(scaling_factor-1.0) > 1E-15:
            if self.echo_level > 2:
                cs_tools.cs_print_info("  Scaling-Factor", scaling_factor)
            interface_data.SetData(scaling_factor*interface_data.GetData()) # setting the scaled data

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "coupling_sequence"        : [],
            "solvers"                  : {},
            "predictors"               : [],
            "coupling_operations"      : {},
            "data_transfer_operators"  : {}
        }""")
        this_defaults.AddMissingParameters(super(CoSimulationCoupledSolver, cls)._GetDefaultSettings())

        return this_defaults

def GetInputDataDefaults():
    return KM.Parameters("""{
        "data"                            : "UNSPECIFIED",
        "from_solver"                     : "UNSPECIFIED",
        "from_solver_data"                : "UNSPECIFIED",
        "data_transfer_operator"          : "UNSPECIFIED",
        "data_transfer_operator_options"  : [],
        "before_data_transfer_operations" : [],
        "after_data_transfer_operations"  : [],
        "interval"                        : [0.0, 1e30],
        "scaling_factor"                  : 1.0
    }""")

def GetOutputDataDefaults():
    return KM.Parameters("""{
        "data"                            : "UNSPECIFIED",
        "to_solver"                       : "UNSPECIFIED",
        "to_solver_data"                  : "UNSPECIFIED",
        "data_transfer_operator"          : "UNSPECIFIED",
        "data_transfer_operator_options"  : [],
        "before_data_transfer_operations" : [],
        "after_data_transfer_operations"  : [],
        "interval"                        : [0.0, 1e30],
        "scaling_factor"                  : 1.0
    }""")
