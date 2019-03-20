from __future__ import print_function, absolute_import, division

# Importing the base class
from  . import co_simulation_base_solver
from collections import OrderedDict

# Other imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import couplingsolverprint, bold

class CoSimulationBaseCouplingSolver(co_simulation_base_solver.CoSimulationBaseSolver):
    def __init__(self, model, cosim_solver_settings):
        super(CoSimulationBaseCouplingSolver, self).__init__(model, cosim_solver_settings)

        # self.solver_names = []
        self.participating_solvers = OrderedDict()

        ### ATTENTION, big flaw, also the participants can be coupled solvers !!!
        import KratosMultiphysics.CoSimulationApplication.solver_interfaces.co_simulation_solver_factory as solver_factory

        for solver_settings in self.cosim_solver_settings["coupling_sequence"]:
            solver_name = solver_settings["name"]
            # if solver_name in self.solver_names:
            #     raise NameError('Solver name "' + solver_name + '" defined twice!')
            # self.solver_names.append(solver_name)
            self.cosim_solver_settings["solvers"][solver_name]["name"] = solver_name # adding the name such that the solver can identify itself
            self.participating_solvers[solver_name] = solver_factory.CreateSolverInterface(
                model, self.cosim_solver_settings["solvers"][solver_name]) # -1 to have solver prints on same lvl

        self.cosim_solver_details = self.__GetSolverCoSimulationDetails(
            self.cosim_solver_settings["coupling_sequence"])

        ### Creating the predictors
        self.predictors_list = cs_tools.CreatePredictors(
            self.cosim_solver_settings["predictors"],
            self.participating_solvers) # TODO set echo-lvl?

        # With this setting the coupling can start later
        self.start_coupling_time = 0.0
        if "start_coupling_time" in self.cosim_solver_settings:
            self.start_coupling_time = self.cosim_solver_settings["start_coupling_time"]
        if self.start_coupling_time > 0.0:
            self.coupling_started = False
        else:
            self.coupling_started = True

    def Initialize(self):
        for solver in self.participating_solvers.values():
            solver.Initialize()
        for solver in self.participating_solvers.values():
            solver.InitializeIO(self.participating_solvers, self.echo_level)
            # we use the Echo_level of the coupling solver, since IO is needed by the coupling
            # and not by the (physics-) solver

        for predictor in self.predictors_list:
            predictor.Initialize()

    def Finalize(self):
        for solver in self.participating_solvers.values():
            solver.Finalize()

        for predictor in self.predictors_list:
            predictor.Finalize()

    def AdvanceInTime(self, current_time):
        self.time = 0.0
        for solver in self.participating_solvers.values():
            self.time = max(self.time, solver.AdvanceInTime(current_time))
        # self.time = self.participating_solvers[self.solver_names[0]].AdvanceInTime(current_time)
        # for solver_name in self.solver_names[1:]:
        #     time_other_solver = self.participating_solvers[solver_name].AdvanceInTime(current_time)
        #     if abs(self.time - time_other_solver) > 1e-12:
        #         raise Exception("Solver time mismatch")

        if not self.coupling_started and self.time > self.start_coupling_time:
            self.coupling_started = True
            if self.echo_level > 0:
                couplingsolverprint(self._Name(), bold("Starting Coupling"))

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

    def FinalizeSolutionStep(self):
        for solver in self.participating_solvers.values():
            solver.FinalizeSolutionStep()
        for predictor in self.predictors_list:
            predictor.FinalizeSolutionStep()

    def OutputSolutionStep(self):
        for solver in self.participating_solvers.values():
            solver.OutputSolutionStep()

    def SolveSolutionStep(self):
        err_msg  = 'Calling "SolveSolutionStep" of the "CoSimulationBaseCouplingSolver"!\n'
        err_msg += 'This function has to be implemented in the derived class!'
        raise Exception(err_msg)

    def _SynchronizeInputData(self, solver, solver_name):
        if self.coupling_started:
            input_data_list = self.cosim_solver_details[solver_name]["input_data_list"]

            if self.time >= self.cosim_solver_details[solver_name]["input_coupling_start_time"]:
                for input_data in input_data_list:
                    from_solver = self.participating_solvers[input_data["from_solver"]]
                    data_name = input_data["data_name"]
                    data_definition = from_solver.GetDataDefinition(data_name)
                    data_settings = { "data_format" : data_definition["data_format"],
                                    "data_name"   : data_name,
                                    "io_settings" : input_data["io_settings"] }
                    solver.ImportCouplingInterfaceData(data_settings, from_solver)

    def _SynchronizeOutputData(self, solver, solver_name):
        if self.coupling_started:
            output_data_list = self.cosim_solver_details[solver_name]["output_data_list"]

            if self.time >= self.cosim_solver_details[solver_name]["output_coupling_start_time"]:
                for output_data in output_data_list:
                    to_solver = self.participating_solvers[output_data["to_solver"]]
                    data_name = output_data["data_name"]
                    data_definition = to_solver.GetDataDefinition(data_name)
                    data_settings = { "data_format" : data_definition["data_format"],
                                    "data_name"   : data_name,
                                    "io_settings" : output_data["io_settings"] }
                    solver.ExportCouplingInterfaceData(data_settings, to_solver)

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


    def __GetSolverCoSimulationDetails(self, co_simulation_solver_settings):
        num_solvers = len(co_simulation_solver_settings)
        solver_cosim_details = {}
        for solver_settings in co_simulation_solver_settings:
            solver_name = solver_settings["name"]
            solver_cosim_details[solver_name] = solver_settings
        return solver_cosim_details
