from __future__ import print_function, absolute_import, division

# Importing the base class
from co_simulation_solvers.co_simulation_base_solver import CoSimulationBaseSolver

# Other imports
from co_simulation_predictors.co_simulation_predictor_factory import CreatePredictor
import co_simulation_tools as cosim_tools
from co_simulation_tools import couplingsolverprint, bold

class CoSimulationBaseCouplingSolver(CoSimulationBaseSolver):
    def __init__(self, cosim_solver_settings, level):
        super(CoSimulationBaseCouplingSolver, self).__init__(cosim_solver_settings, level)

        self.solver_names = []
        self.solvers = {}

        import python_solvers_wrapper_co_simulation as solvers_wrapper

        for solver_settings in self.cosim_solver_settings["coupling_loop"]:
            solver_name = solver_settings["name"]
            if solver_name in self.solver_names:
                raise NameError('Solver name "' + solver_name + '" defined twice!')
            self.solver_names.append(solver_name)
            self.cosim_solver_settings["solvers"][solver_name]["name"] = solver_name # adding the name such that the solver can identify itself
            self.solvers[solver_name] = solvers_wrapper.CreateSolver(
                self.cosim_solver_settings["solvers"][solver_name], self.lvl-1) # -1 to have solver prints on same lvl

        self.cosim_solver_details = cosim_tools.GetSolverCoSimulationDetails(
            self.cosim_solver_settings["coupling_loop"])

        self.predictor = None
        if "predictor_settings" in self.cosim_solver_settings:
            self.predictor = CreatePredictor(self.cosim_solver_settings["predictor_settings"],
                                             self.solvers, self.lvl)
            self.predictor.SetEchoLevel(self.echo_level)

        # With this setting the coupling can start later
        self.start_coupling_time = 0.0
        if "start_coupling_time" in self.cosim_solver_settings:
            self.start_coupling_time = self.cosim_solver_settings["start_coupling_time"]
        if self.start_coupling_time > 0.0:
            self.coupling_started = False
        else:
            self.coupling_started = True

    def Initialize(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].Initialize()
        for solver_name in self.solver_names:
            self.solvers[solver_name].InitializeIO(self.solvers, self.echo_level)
            # we use the Echo_level of the coupling solver, since IO is needed by the coupling
            # and not by the (physics-) solver

        if self.predictor is not None:
            self.predictor.Initialize()

    def Finalize(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].Finalize()

        if self.predictor is not None:
            self.predictor.Finalize()

    def AdvanceInTime(self, current_time):
        self.time = self.solvers[self.solver_names[0]].AdvanceInTime(current_time)
        for solver_name in self.solver_names[1:]:
            time_other_solver = self.solvers[solver_name].AdvanceInTime(current_time)
            if abs(self.time - time_other_solver) > 1e-12:
                raise Exception("Solver time mismatch")

        if not self.coupling_started and self.time > self.start_coupling_time:
            self.coupling_started = True
            if self.echo_level > 0:
                couplingsolverprint(self.lvl, self._Name(), bold("Starting Coupling"))

        # if a predictor is used then the delta_time is set
        # this is needed by some predictors
        if self.predictor is not None:
            delta_time = self.time - current_time
            self.predictor.SetDeltaTime(delta_time)

        return self.time

    def Predict(self):
        if self.predictor is not None:
            self.predictor.Predict()
        for solver_name in self.solver_names:
            self.solvers[solver_name].Predict()

    def InitializeSolutionStep(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].InitializeSolutionStep()

        if self.predictor is not None:
            self.predictor.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].FinalizeSolutionStep()
        if self.predictor is not None:
            self.predictor.FinalizeSolutionStep()

    def OutputSolutionStep(self):
        for solver_name in self.solver_names:
            self.solvers[solver_name].OutputSolutionStep()

    def SolveSolutionStep(self):
        err_msg  = 'Calling "SolveSolutionStep" of the "CoSimulationBaseCouplingSolver"!\n'
        err_msg += 'This function has to be implemented in the derived class!'
        raise Exception(err_msg)

    def _SynchronizeInputData(self, solver, solver_name):
        if self.coupling_started:
            input_data_list = self.cosim_solver_details[solver_name]["input_data_list"]
            for input_data in input_data_list:
                from_solver = self.solvers[input_data["from_solver"]]
                data_name = input_data["data_name"]
                data_definition = from_solver.GetDataDefinition(data_name)
                data_settings = { "data_format" : data_definition["data_format"],
                                  "data_name"   : data_name,
                                  "io_settings" : input_data["io_settings"] }
                solver.ImportData(data_settings, from_solver)

    def _SynchronizeOutputData(self, solver, solver_name):
        if self.coupling_started:
            output_data_list = self.cosim_solver_details[solver_name]["output_data_list"]
            for output_data in output_data_list:
                to_solver = self.solvers[output_data["to_solver"]]
                data_name = output_data["data_name"]
                data_definition = to_solver.GetDataDefinition(data_name)
                data_settings = { "data_format" : data_definition["data_format"],
                                  "data_name"   : data_name,
                                  "io_settings" : output_data["io_settings"] }
                solver.ExportData(data_settings, to_solver)

    def PrintInfo(self):
        super(CoSimulationBaseCouplingSolver, self).PrintInfo()

        couplingsolverprint(self.lvl, self._Name(), "Has the following participants:")

        for solver_name in self.solver_names:
            self.solvers[solver_name].PrintInfo()

        if self.predictor is not None:
            couplingsolverprint(self.lvl, self._Name(), "Uses a Predictor:")
            self.predictor.PrintInfo()

    def Check(self):
        super(CoSimulationBaseCouplingSolver, self).Check()

        for solver_name in self.solver_names:
            self.solvers[solver_name].Check()

        if self.predictor is not None:
            self.predictor.Check()

    def IsDistributed(self):
        return True
