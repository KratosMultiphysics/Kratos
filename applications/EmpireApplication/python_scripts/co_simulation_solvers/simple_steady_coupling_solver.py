# Importing the base class
from co_simulation_solvers.gauss_seidel_strong_coupling_solver import GaussSeidelStrongCouplingSolver

def CreateSolver(cosim_solver_settings, level):
    return SimpleSteadyCouplingSolver(cosim_solver_settings, level)

class SimpleSteadyCouplingSolver(GaussSeidelStrongCouplingSolver):

    def _SynchronizeInputData(self, solver, solver_name):

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
        output_data_list = self.cosim_solver_details[solver_name]["output_data_list"]

        for output_data in output_data_list:
            to_solver = self.solvers[output_data["to_solver"]]
            data_name = output_data["data_name"]
            data_definition = to_solver.GetDataDefinition(data_name)
            data_settings = { "data_format" : data_definition["data_format"],
                            "data_name"   : data_name,
                            "io_settings" : output_data["io_settings"] }
            solver.ExportData(data_settings, to_solver)