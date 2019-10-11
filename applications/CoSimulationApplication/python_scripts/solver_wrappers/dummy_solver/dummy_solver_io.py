from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_io import CoSimulationIO

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication as KratosCoSim
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

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

        self.io = KratosCoSim.CoSimIO(solver_name, ParametersToStringDict(self.settings.items()))

        self.io.Connect()

    def Finalize(self):
        self.io.Disconnect()

    # def ImportCouplingInterface(self, interface_config):
    #     model_part_name = interface_config["model_part_name"]
    #     comm_name = interface_config["comm_name"]

    #     if not self.model.HasModelPart(model_part_name):
    #         main_model_part_name, *sub_model_part_names = model_part_name.split(".")
    #         cs_tools.RecursiveCreateModelParts(self.model[main_model_part_name], ".".join(sub_model_part_names))

    #     model_part = self.model[model_part_name]
    #     KratosCoSim.EMPIRE_API.EMPIRE_API_recvMesh(model_part, comm_name)

    # def ExportCouplingInterface(self, interface_config):
    #     model_part_name = interface_config["model_part_name"]
    #     comm_name = interface_config["comm_name"]
    #     KratosCoSim.EMPIRE_API.EMPIRE_API_sendMesh(self.model[model_part_name], comm_name)

    # def ImportData(self, data_config):
    #     data_type = data_config["type"]
    #     if data_type == "coupling_interface_data":
    #         interface_data = data_config["interface_data"]
    #         KratosCoSim.EMPIRE_API.EMPIRE_API_recvDataField(interface_data.GetModelPart(), self.solver_name+"_"+interface_data.name, interface_data.variable)
    #     else:
    #         raise NotImplementedError('Importing interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def ExportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            self.ExportData(interface_data.GetData()) # TODO come up with sth better here, this copies every time!
        elif data_type == "control_signal":
            control_signal_key = cs.Tools.control_signal_map[data_config["signal"]]
            self.io.ExportData(control_signal_key)
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
