from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_io import CoSimulationBaseIO

# Other imports
import numpy as np
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(solvers, solver_name):
    return SDoFIO(solvers, solver_name)

class SDoFIO(CoSimulationBaseIO):

    def ImportCouplingInterfaceData(self, data_settings, from_client):
        data_name = data_settings["data_name"]
        io_settings = data_settings["io_settings"]
        data_array = np.array([])
        cs_tools.ImportArrayFromSolver(from_client, data_name, data_array)

        sdof_solver = self.solvers[self.solver_name]
        sdof_data_settings = sdof_solver.GetInterfaceData(data_settings["data_name"])

        value = sum(data_array)

        if io_settings.Has("io_options"):
            if io_settings["io_options"].Has("swap_sign"):
                value *= -1.0

        data_identifier = sdof_data_settings["data_identifier"].GetString()
        sdof_solver.SetData(data_identifier, value)


    def ExportCouplingInterfaceData(self, data_settings, to_client):
        sdof_solver = self.solvers[self.solver_name]
        sdof_data_settings = sdof_solver.GetInterfaceData(data_settings["data_name"])
        out_data_settings = {"data_format" : "scalar_value"}
        out_data_settings["data_name"] = data_settings["data_name"]
        # sdof_data_settings["data_name"].SetString(data_settings["data_name"])# TODO maybe it is not existing before ...? => check
        # sdof_data_settings.AddEmptyValue("data_name").SetString(data_settings["data_name"])# TODO maybe it is not existing before ...? => check

        if not sdof_data_settings["data_format"].GetString() == "scalar_value":
            raise Exception('SDoFIO can only handle scalar values')
        sdof_solver = self.solvers[self.solver_name]

        data_identifier = sdof_data_settings["data_identifier"].GetString()

        x = sdof_solver.GetData(data_identifier)

        out_data_settings["scalar_value"] = x
        # sdof_data_settings.AddEmptyValue("scalar_value").SetDouble(x)
        # print(sdof_data_settings)
        # TODO convert to dict
        to_client.ImportCouplingInterfaceData(data_settings, to_client)
