from __future__ import print_function, absolute_import, division

# Importing the base class
from  . import co_simulation_base_solver

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import couplingsolverprint, bold

class CoSimulationBaseCouplingSolver(co_simulation_base_solver.CoSimulationBaseSolver):
    def __init__(self, model, cosim_solver_settings, solver_name):
        super(CoSimulationBaseCouplingSolver, self).__init__(model, cosim_solver_settings, solver_name)

        self.participating_solvers = self.__CreateSolvers()
        self.coupling_sequence = self.__GetSolverCoSimulationDetails()

        ### Creating the predictors
        self.predictors_list = cs_tools.CreatePredictors(
            self.settings["predictors"],
            self.participating_solvers,
            self.echo_level)

        ### Creating the coupling operations
        self.coupling_operations_dict = cs_tools.CreateCouplingOperations(
            self.settings["coupling_operations"],
            self.participating_solvers,
            self.echo_level)

        self.coupling_operations_list = [d for d in self.coupling_operations_dict.values()]

        ### Creating the data transfer operators
        self.data_transfer_operators_dict = cs_tools.CreateDataTransferOperators(
            self.settings["data_transfer_operators"],
            self.echo_level)

    def Initialize(self):
        for solver in self.participating_solvers.values():
            solver.Initialize()
        for solver in self.participating_solvers.values():
            solver.InitializeIO(self.participating_solvers, self.echo_level)
            # we use the Echo_level of the coupling solver, since IO is needed by the coupling
            # and not by the (physics-) solver

        for predictor in self.predictors_list:
            predictor.Initialize()

        for coupling_op in self.coupling_operations_list:
            coupling_op.Initialize()

    def Finalize(self):
        for solver in self.participating_solvers.values():
            solver.Finalize()

        for predictor in self.predictors_list:
            predictor.Finalize()

        for coupling_op in self.coupling_operations_list:
            coupling_op.Finalize()

    def AdvanceInTime(self, current_time):
        self.time = 0.0
        for solver in self.participating_solvers.values():
            self.time = max(self.time, solver.AdvanceInTime(current_time))

        return self.time

    def Predict(self):
        for predictor in self.predictors_list:
            predictor.Predict()

        for solver in self.participating_solvers.values():
            solver.Predict()

    def InitializeSolutionStep(self):
        for solver in self.participating_solvers.values():
            solver.InitializeSolutionStep()

        for predictor in self.predictors_list:
            predictor.InitializeSolutionStep()

        for coupling_op in self.coupling_operations_list:
            coupling_op.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        for solver in self.participating_solvers.values():
            solver.FinalizeSolutionStep()

        for predictor in self.predictors_list:
            predictor.FinalizeSolutionStep()

        for coupling_op in self.coupling_operations_list:
            coupling_op.FinalizeSolutionStep()

    def OutputSolutionStep(self):
        for solver in self.participating_solvers.values():
            solver.OutputSolutionStep()

    def SolveSolutionStep(self):
        err_msg  = 'Calling "SolveSolutionStep" of the "CoSimulationBaseCouplingSolver"!\n'
        err_msg += 'This function has to be implemented in the derived class!'
        raise Exception(err_msg)

    def _SynchronizeInputData(self, solver_name):
        to_solver = self.participating_solvers[solver_name]
        input_data_list = self.coupling_sequence[solver_name]["input_data_list"]

        # TODO check if in interval
        for i in range(input_data_list.size()):
            i_input_data = input_data_list[i]

            # from solver
            from_solver = self.participating_solvers[i_input_data["from_solver"].GetString()]
            from_solver_data = from_solver.GetInterfaceData(i_input_data["from_solver_data"].GetString())

            # to solver
            to_solver_data = to_solver.GetInterfaceData(i_input_data["data"].GetString())

            # perform the data transfer
            self.__ExecuteCouplingOperations(i_input_data["before_data_transfer_operations"])

            data_transfer_operator_name = i_input_data["data_transfer_operator"].GetString()
            # TODO check the order of solvers!
            self.__GetDataTransferOperator(data_transfer_operator_name).TransferData(from_solver, to_solver, i_input_data["data_transfer_operator_options"])

            self.__ExecuteCouplingOperations(i_input_data["after_data_transfer_operations"])
            if i_input_data.Has("mapper_settings"):
                to_solver_data.origin_data = from_solver_data
                to_solver_data.mapper_settings = i_input_data["mapper_settings"]

            to_solver.ImportCouplingInterfaceData(to_solver_data, from_solver)


    def _SynchronizeOutputData(self, solver_name):
        from_solver = self.participating_solvers[solver_name]
        output_data_list = self.coupling_sequence[solver_name]["output_data_list"]

        for i in range(output_data_list.size()):
            i_output_data = output_data_list[i]

            # from solver
            from_solver_data = from_solver.GetInterfaceData(i_output_data["data"].GetString())

            # to solver
            to_solver = self.participating_solvers[i_output_data["to_solver"].GetString()]
            to_solver_data = to_solver.GetInterfaceData(i_output_data["to_solver_data"].GetString())

            # perform the data transfer
            self.__ExecuteCouplingOperations(i_output_data["before_data_transfer_operations"])

            data_transfer_operator_name = i_output_data["data_transfer_operator"].GetString()
            # TODO check the order of solvers!
            self.__GetDataTransferOperator(data_transfer_operator_name).TransferData(from_solver, to_solver, i_output_data["data_transfer_operator_options"])

            self.__ExecuteCouplingOperations(i_output_data["after_data_transfer_operations"])

            if i_output_data.Has("mapper_settings"):
                from_solver_data.destination_data = to_solver_data
                from_solver_data.mapper_settings = i_output_data["mapper_settings"]

            # export the data
            from_solver.ExportCouplingInterfaceData(from_solver_data, to_solver) # TODO I am quite sure that the args have to be changed



    def __GetDataTransferOperator(self, data_transfer_operator_name):
        try:
            return self.data_transfer_operators_dict[data_transfer_operator_name]
        except KeyError:
            raise NameError('The data-transfer-operator "{}" does not exist!'.format(data_transfer_operator_name))


    def __ExecuteCouplingOperations(self, settings):
        for i in range(settings.size()):
            name = settings[i].GetString()
            self.coupling_operations_dict[coupling_operation_name].Execute()

    def PrintInfo(self):
        super(CoSimulationBaseCouplingSolver, self).PrintInfo()

        couplingsolverprint(self._Name(), "Has the following participants:")

        for solver in self.participating_solvers.values():
            solver.PrintInfo()

        # if self.predictor is not None:
        #     couplingsolverprint(self._Name(), "Uses a Predictor:")
        #     self.predictor.PrintInfo()

    def Check(self):
        super(CoSimulationBaseCouplingSolver, self).Check()

        for solver in self.participating_solvers.values():
            solver.Check()

        # if self.predictor is not None:
        #     self.predictor.Check()

    def IsDistributed(self):
        return True

    def _GetIOName(self):
        # the coupled-solvers always use the kratos-format, since there are no "external" coupled solvers
        return "kratos"

    def __CreateSolvers(self):
        ### ATTENTION, big flaw, also the participants can be coupled solvers !!!
        import KratosMultiphysics.CoSimulationApplication.solver_interfaces.co_simulation_solver_factory as solver_factory
        from collections import OrderedDict
        # first create all solvers
        solvers = {}
        for solver_name, solver_settings in self.settings["solvers"].items():
            solvers[solver_name] = solver_factory.CreateSolverInterface(
                self.model, solver_settings, solver_name)

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
                data_list[i_data_list].ValidateAndAssignDefaults(defaults)

        solver_cosim_details = {}
        for i_solver_settings in range(self.settings["coupling_sequence"].size()):
            solver_settings = self.settings["coupling_sequence"][i_solver_settings]
            solver_name = solver_settings["name"].GetString()
            solver_cosim_details[solver_name] = solver_settings

            ValidateAndAssignDefaultsDataList(solver_settings["input_data_list"], GetInputDataDefaults())
            ValidateAndAssignDefaultsDataList(solver_settings["output_data_list"], GetOutputDataDefaults())

        return solver_cosim_details

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = cs_tools.cs_data_structure.Parameters("""{
            "coupling_sequence"        : [],
            "solvers"                  : {},
            "predictors"               : [],
            "coupling_operations"      : {},
            "data_transfer_operators"  : {}
        }""")
        this_defaults.AddMissingParameters(super(CoSimulationBaseCouplingSolver, cls)._GetDefaultSettings())

        return this_defaults

def GetInputDataDefaults():
    return cs_tools.cs_data_structure.Parameters("""{
        "data"                            : "UNSPECIFIED",
        "from_solver"                     : "UNSPECIFIED",
        "from_solver_data"                : "UNSPECIFIED",
        "data_transfer_operator"          : "UNSPECIFIED",
        "data_transfer_operator_options"  : [],
        "before_data_transfer_operations" : [],
        "after_data_transfer_operations"  : [],
        "interval"                        : [],
        "mapper_settings" : {}
    }""")

def GetOutputDataDefaults():
    return cs_tools.cs_data_structure.Parameters("""{
        "data"                            : "UNSPECIFIED",
        "to_solver"                       : "UNSPECIFIED",
        "to_solver_data"                  : "UNSPECIFIED",
        "data_transfer_operator"          : "UNSPECIFIED",
        "data_transfer_operator_options"  : [],
        "before_data_transfer_operations" : [],
        "after_data_transfer_operations"  : [],
        "interval"                        : [],
        "mapper_settings" : {}
    }""")
