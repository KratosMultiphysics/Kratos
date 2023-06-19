# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.factories.solver_wrapper_factory as solver_wrapper_factory
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.factories.helpers as factories_helper
import KratosMultiphysics.CoSimulationApplication.colors as colors
from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import BaseCouplingInterfaceData

# Other imports
from collections import OrderedDict

class UndefinedSolver:
    def __init__(self, name, settings):
        self.name = name
        self.settings = settings

    def Initialize(self):
        if self.settings.Has("data"):
            self.data_dict = {data_name : BaseCouplingInterfaceData(data_config, data_name, self.name) for (data_name, data_config) in self.settings["data"].items()}
        else:
            self.data_dict = {}

    def IsDefinedOnThisRank(self):
        return False

    def GetInterfaceData(self, data_name):
        try:
            return self.data_dict[data_name]
        except KeyError:
            raise Exception('Requested data field "{}" does not exist for solver "{}"'.format(data_name, self.name))

    def AdvanceInTime(*args): return 0.0

    def __getattr__(self, attr):
        return lambda *args : None

class CoSimulationCoupledSolver(CoSimulationSolverWrapper):
    """Baseclass for the coupled solvers used for CoSimulation
    Performs basic operations that are common among coupled solvers:
    - holds Predictors
    - holds DataTransferOperators
    - holds CouplingOperations
    - initialization of IOs of solvers
    - Synchronization of Input and Output
    - Handles the coupling sequence
    """
    def __init__(self, settings, models, solver_name):
        # perform some initial checks
        if not settings.Has("coupling_sequence"):
            err_msg  = 'No "coupling_sequence" was specified for coupled solver\n'
            err_msg += '"{}" of type "{}"'.format(solver_name, self._ClassName())
            raise Exception(err_msg)

        if settings["coupling_sequence"].size() == 0:
            err_msg  = '"coupling_sequence" is empty for coupled solver\n'
            err_msg += '"{}" of type "{}"'.format(solver_name, self._ClassName())
            raise Exception(err_msg)

        if not settings.Has("solvers"):
            err_msg  = 'No "solvers" are specified for coupled solver\n'
            err_msg += '"{}" of type "{}"'.format(solver_name, self._ClassName())
            raise Exception(err_msg)

        if len(settings["solvers"].keys()) == 0:
            err_msg  = '"solvers" is empty for coupled solver\n'
            err_msg += '"{}" of type "{}"'.format(solver_name, self._ClassName())
            raise Exception(err_msg)

        if not isinstance(models, dict) and not models is None:
            err_msg  = 'A coupled solver can either be passed a dict of Models\n'
            err_msg += 'or None, got object of type "{}"'.format(type(models))
            raise Exception(err_msg)

        super().__init__(settings, None, solver_name)

        self.process_info = KM.ProcessInfo()

        # TODO initialize this in a restart
        self.process_info[KM.STEP] = 0
        self.process_info[KM.TIME] = 0.0
        self.process_info[KM.IS_RESTARTED] = False

        self.solver_wrappers = self.__CreateSolverWrappers(models)

        # overwriting the Model created in the BaseClass
        # CoupledSolvers only forward calls to its solvers
        # this is done with the ModelAccessor
        self.model = ModelAccessor(self.solver_wrappers)

        self.coupling_sequence = self.__GetSolverCoSimulationDetails()

        for solver in self.solver_wrappers.values():
            solver.CreateIO(self.echo_level)
            # using the echo_level of the coupled solver, since IO is needed by the coupling

    def _GetSolver(self, solver_name):
        solver_name, *sub_solver_names = solver_name.split(".")
        solver = self.solver_wrappers[solver_name]
        if len(sub_solver_names) > 0:
            return solver._GetSolver(".".join(sub_solver_names))
        else:
            return solver

    def Initialize(self):
        for solver in self.solver_wrappers.values():
            solver.Initialize()

        super().Initialize()

        ### Creating the predictors
        self.predictors_list = factories_helper.CreatePredictors(
            self.settings["predictors"],
            self.solver_wrappers,
            self.echo_level)

        ### Creating the coupling operations
        self.coupling_operations_dict = factories_helper.CreateCouplingOperations(
            self.settings["coupling_operations"],
            self.solver_wrappers,
            self.process_info,
            self.data_communicator,
            self.echo_level)

        ### Creating the data transfer operators
        self.data_transfer_operators_dict = factories_helper.CreateDataTransferOperators(
            self.settings["data_transfer_operators"],
            self.data_communicator,
            self.echo_level)

        for predictor in self.predictors_list:
            predictor.Initialize()

        for coupling_operation in self.coupling_operations_dict.values():
            coupling_operation.Initialize()

    def Finalize(self):
        super().Finalize()

        for solver in self.solver_wrappers.values():
            solver.Finalize()

        for predictor in self.predictors_list:
            predictor.Finalize()

        for coupling_operation in self.coupling_operations_dict.values():
            coupling_operation.Finalize()

    def AdvanceInTime(self, current_time):
        # not all solvers provide time (e.g. external solvers or steady solvers)
        # hence we have to check first if they return time (i.e. time != 0.0)
        # and then if the times are matching, since currently no interpolation in time is possible

        self.time = 0.0
        for solver in self.solver_wrappers.values():
            # TODO maybe do a check to make sure all ranks have the same time?
            solver_time = self.data_communicator.MaxAll(solver.AdvanceInTime(current_time))
            if solver_time != 0.0: # solver provides time
                if self.time == 0.0: # first time a solver returns a time different from 0.0
                    self.time = solver_time
                elif abs(self.time - solver_time) > 1e-12:
                        raise Exception("Solver time mismatch")

        self.process_info[KM.TIME] = self.time
        self.process_info[KM.STEP] += 1

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
        data_list = self.coupling_sequence[solver_name]["input_data_list"]
        if self.echo_level > 2:
            cs_tools.cs_print_info(self._ClassName(), 'Start Synchronizing Input for solver "{}"'.format(colors.blue(solver_name)))

        for i in range(data_list.size()):
            i_data = data_list[i]

            to_data_name = i_data["data"].GetString()
            from_solver_name = i_data["from_solver"].GetString()
            from_solver_data_name = i_data["from_solver_data"].GetString()

            if self.echo_level > 2:
                cs_tools.cs_print_info("  Data", '"{}" | from solver: "{}": "{}"'.format(colors.magenta(to_data_name), colors.blue(from_solver_name), colors.magenta(from_solver_data_name)))

            # from solver
            from_solver_data = self.__GetInterfaceDataFromSolver(from_solver_name, from_solver_data_name)

            # to solver
            to_solver_data = self.__GetInterfaceDataFromSolver(solver_name, to_data_name)

            self.__SynchronizeData(i_data, from_solver_data, to_solver_data)

        if self.echo_level > 2:
            cs_tools.cs_print_info(self._ClassName(), 'End Synchronizing Input for solver "{}"'.format(colors.blue(solver_name)))

    def _SynchronizeOutputData(self, solver_name):
        data_list = self.coupling_sequence[solver_name]["output_data_list"]
        if self.echo_level > 2:
            cs_tools.cs_print_info(self._ClassName(), 'Start Synchronizing Output for solver "{}"'.format(colors.blue(solver_name)))

        for i in range(data_list.size()):
            i_data = data_list[i]

            from_data_name = i_data["data"].GetString()
            to_solver_name = i_data["to_solver"].GetString()

            to_solver_data_name = i_data["to_solver_data"].GetString()

            if self.echo_level > 2:
                cs_tools.cs_print_info("  Data", '"{}" | to solver: "{}": "{}"'.format(colors.magenta(from_data_name), colors.blue(to_solver_name), colors.magenta(to_solver_data_name)))

            # from solver
            from_solver_data = self.__GetInterfaceDataFromSolver(solver_name, from_data_name)

            # to solver
            to_solver_data = self.__GetInterfaceDataFromSolver(to_solver_name, to_solver_data_name)

            self.__SynchronizeData(i_data, from_solver_data, to_solver_data)

        if self.echo_level > 2:
            cs_tools.cs_print_info(self._ClassName(), 'End Synchronizing Output for solver "{}"'.format(colors.blue(solver_name)))

    def __SynchronizeData(self, i_data, from_solver_data, to_solver_data):
        # Check if data-exchange is specified for current time
        if not KM.IntervalUtility(i_data).IsInInterval(self.time):
            if self.echo_level > 2:
                cs_tools.cs_print_info("  Skipped", 'not in interval')
            return

        # Perform the data transfer
        self.__ExecuteCouplingOperations(i_data["before_data_transfer_operations"])

        data_transfer_operator_name = i_data["data_transfer_operator"].GetString()
        # TODO check the order of solvers!
        self.__GetDataTransferOperator(data_transfer_operator_name).TransferData(from_solver_data, to_solver_data, i_data["data_transfer_operator_options"])

        self.__ExecuteCouplingOperations(i_data["after_data_transfer_operations"])

    def __GetInterfaceDataFromSolver(self, solver_name, interface_data_name):
        solver = self.solver_wrappers[solver_name]
        return solver.GetInterfaceData(interface_data_name)

    def __GetDataTransferOperator(self, data_transfer_operator_name):
        try:
            return self.data_transfer_operators_dict[data_transfer_operator_name]
        except KeyError:
            raise NameError('The data-transfer-operator "{}" does not exist!'.format(data_transfer_operator_name))

    def __ExecuteCouplingOperations(self, settings):
        for coupling_operation_name in settings.GetStringArray():
            self.coupling_operations_dict[coupling_operation_name].Execute()

    def PrintInfo(self):
        super().PrintInfo()

        cs_tools.cs_print_info(self._ClassName(), "Has the following components:")
        for solver in self.solver_wrappers.values():
            solver.PrintInfo()

        for predictor in self.predictors_list:
            predictor.PrintInfo()

        for coupling_operation in self.coupling_operations_dict.values():
            coupling_operation.PrintInfo()

    def Check(self):
        # TODO check that there is no self-communication with the same data!
        # self-communication is allowed within a solver, but not on the same data
        super().Check()

        for solver in self.solver_wrappers.values():
            solver.Check()

        for predictor in self.predictors_list:
            predictor.Check()

        for coupling_operation in self.coupling_operations_dict.values():
            coupling_operation.Check()

    def __CreateSolverWrappers(self, models):
        # first create all solvers
        solvers = {}
        for solver_name, solver_settings in self.settings["solvers"].items():
            if models == None:
                solver_model = None
            else:
                solver_model = models.get(solver_name) # returns None if "solver_name" is not in models
            solvers[solver_name] = solver_wrapper_factory.CreateSolverWrapper(solver_settings, solver_model, solver_name)

        # then order them according to the coupling-loop
        solvers_map = OrderedDict()
        for i_solver_settings in range(self.settings["coupling_sequence"].size()):
            solver_settings = self.settings["coupling_sequence"][i_solver_settings]
            solver_name = solver_settings["name"].GetString()
            solver = solvers[solver_name]
            if solver.IsDefinedOnThisRank():
                solvers_map[solver_name] = solvers[solver_name]
            else:
                solvers_map[solver_name] = UndefinedSolver(solver_name, self.settings["solvers"][solver_name])

        for solver_name in self.settings["solvers"].keys():
            if solver_name not in solvers_map:
                err_msg  = 'Solver "{}" of type "{}"\n'.format(solver_name, solvers[solver_name]._ClassName())
                err_msg += 'is specified in the "solvers" of coupled solver\n'
                err_msg += '"{}" of type "{}"\n'.format(self.name, self._ClassName())
                err_msg += 'but not used in the "coupling_sequence"!'
                raise Exception(err_msg)

        if models != None:
            for solver_name in models.keys():
                if solver_name not in solvers_map:
                    raise Exception('A Model was given for solver "{}" but this solver does not exist!'.format(solver_name))

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

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "coupling_sequence"        : [],
            "solvers"                  : {},
            "predictors"               : [],
            "coupling_operations"      : {},
            "data_transfer_operators"  : {}
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())

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
        "interval"                        : [0.0, 1e30]
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
        "interval"                        : [0.0, 1e30]
    }""")


class ModelAccessor(object):
    """Intermediate class for redirecting the access to the Models
    to the solvers of the CoupledSolver
    """
    def __init__(self, solver_wrappers):
        self.solver_wrappers = solver_wrappers

    def __getitem__(self, key):
        splitted_key = key.split('.')

        if key == "":
            raise Exception("No solver_name was specified!")
        elif key.count('.') == 0:
            # if only the solver name was given then return the Model itself
            return self.solver_wrappers[splitted_key[0]].model

        solver_name, model_part_name = key.split('.', 1)
        # note that model_part_name can still include solver-names in a multicoupling scenario

        return self.solver_wrappers[solver_name].model[model_part_name]
