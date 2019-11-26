from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_io import CoSimulationIO

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.CoSimIO as CoSimIO

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

        KratosCoSimCoSimIO.Connect(self.solver_name, ParametersToStringDict(self.settings), True)

    def Finalize(self):
        KratosCoSimCoSimIO.Disconnect(self.solver_name)

    def ImportCouplingInterface(self, interface_config):
        CoSimIO.ImportMesh(self.solver_name, self.model[interface_config["model_part_name"]]) # TODO this can also be geometry at some point

    def ImportCouplingInterface(self, interface_config):
        CoSimIO.ExportMesh(self.solver_name, self.model[interface_config["model_part_name"]]) # TODO this can also be geometry at some point

    def ImportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            CoSimIO.ImportData(self.solver_name, interface_data.name, interface_data.GetModelPart(), interface_data.variable, GetDataLocation(interface_data.location))

        elif data_type == "time":
            time_list = [0.0]
            KratosCoSim.CoSimIO.ImportData("time_to_co_sim", time_list)
            data_config["time"] = time_list[0]
        else:
            raise NotImplementedError('Exporting interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def ExportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            CoSimIO.ExportData(self.solver_name, interface_data.name, interface_data.GetModelPart(), interface_data.variable, GetDataLocation(interface_data.location))

        elif data_type == "control_signal":
            control_signal_key = data_config["signal"]
            CoSimIO.SendControlSignal(self.solver_name, data_config["identifier"], control_signal_key)

        elif data_type == "time":
            current_time = data_config["time"]
            CoSimIO.ExportData(self.solver_name, "time_from_co_sim", [current_time])

        elif data_type == "convergence_signal":
            if data_config["is_converged"]:
                control_signal_key = CoSimIO.ControlSignal.ConvergenceAchieved
            else:
                control_signal_key = CoSimIO.ControlSignal.Dummy
            CoSimIO.SendControlSignal(self.solver_name, "", control_signal_key)
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
        "node_historical"     : CoSimIO.DataLocation.NodeHistorical,
        "node_non_historical" : CoSimIO.DataLocation.NodeNonHistorical,
        "element"             : CoSimIO.DataLocation.Element,
        "condition"           : CoSimIO.DataLocation.Condition,
        "model_part"          : CoSimIO.DataLocation.ModelPart
    }
    return location_map[location_str]
