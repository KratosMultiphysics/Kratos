from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_io import CoSimulationIO

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication as KratosCoSim

# Other imports
import os

def Create(model, settings, solver_name):
    return DummySolverIO(model, settings, solver_name)

class DummySolverIO(CoSimulationIO):
    """This class is used if a Solver directly uses Kratos as a data-structure
    e.g. Kratos itself or simple-solvers written in Python
    """
    def __init__(self, settings, model, solver_name):
        super(DummySolverIO, self).__init__(settings, model, solver_name)

        self.io = KratosCoSim.CoSimIO(solver_name, ParametersToStringDict(self.settings), True)

        self.io.Connect()

    def Finalize(self):
        self.io.Disconnect()

    def ImportCouplingInterface(self, interface_config):
        self.io.ImportMesh(self.model[interface_config["model_part_name"]]) # TODO this can also be geometry at some point

    def ImportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            self.io.ImportData(interface_data.GetModelPart(), interface_data.variable, GetDataLocation(interface_data.location), interface_data.name)
        else:
            raise NotImplementedError('Exporting interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def ExportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            self.io.ExportData(interface_data.GetModelPart(), interface_data.variable, GetDataLocation(interface_data.location), interface_data.name)

        elif data_type == "control_signal":
            control_signal_key = data_config["signal"]
            self.io.SendControlSignal(control_signal_key, data_config["identifier"])

        elif data_type == "convergence_signal":
            control_signal_key = 0
            if data_config["is_converged"]:
                control_signal_key = 51
            self.io.SendControlSignal(control_signal_key, "dummy")
        else:
            raise NotImplementedError('Exporting interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def PrintInfo(self):
        print("This is the EMPIRE-IO")

    def Check(self):
        pass


def ParametersToStringDict(param):
    string_dict = {}

    for k,v in param.items():
        if v.IsInt():
            v = v.GetInt()
        elif v.IsBool():
            v = int(v.GetBool())
        elif v.IsDouble():
            v = v.GetDouble()
        elif v.IsString():
            v = v.GetString()
        else:
            raise Exception("Only int, double, bool and string are allowed!")

        string_dict[k] = str(v)

    return string_dict

def GetDataLocation(location_str):
    location_map = {
        "node_historical" : KratosCoSim.DataLocation.NodeHistorical,
        "node_non_historical" : KratosCoSim.DataLocation.NodeNonHistorical,
        "element" : KratosCoSim.DataLocation.Element,
        "condition" : KratosCoSim.DataLocation.Condition,
        "model_part" : KratosCoSim.DataLocation.ModelPart
    }
    return location_map[location_str]
