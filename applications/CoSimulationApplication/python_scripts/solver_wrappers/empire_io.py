from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics import kratos_utilities

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_io import CoSimulationIO

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication as KratosCoSim
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

# Other imports
import os

def Create(settings, model, solver_name):
    return EmpireIO(settings, model, solver_name)

class EmpireIO(CoSimulationIO):
    """IO for the legacy EMPIRE_API
    """
    def __init__(self, settings, model, solver_name):
        super(EmpireIO, self).__init__(settings, model, solver_name)
        KratosCoSim.EMPIRE_API.EMPIRE_API_Connect(self.settings["api_configuration_file_name"].GetString())

        # delete and recreate communication folder to avoid leftover files
        self.communication_folder = self.settings["communication_folder"].GetString()
        kratos_utilities.DeleteDirectoryIfExisting(self.communication_folder)
        os.mkdir(self.communication_folder)

        self.double_vector_var_map = {KM.KratosGlobals.GetVariable(k) : KM.KratosGlobals.GetVariable(v.GetString()), in (k,v) in self.settings["communicate_as_double_vector"].items()}

    def Finalize(self):
        kratos_utilities.DeleteDirectoryIfExisting(self.communication_folder)

    def ImportCouplingInterface(self, interface_config):
        model_part_name = interface_config["model_part_name"]
        comm_name = interface_config["comm_name"]

        if not self.model.HasModelPart(model_part_name):
            main_model_part_name, *sub_model_part_names = model_part_name.split(".")
            cs_tools.RecursiveCreateModelParts(self.model[main_model_part_name], ".".join(sub_model_part_names))

        model_part = self.model[model_part_name]
        KratosCoSim.EMPIRE_API.EMPIRE_API_recvMesh(model_part, comm_name)

    def ExportCouplingInterface(self, interface_config):
        model_part_name = interface_config["model_part_name"]
        comm_name = interface_config["comm_name"]
        KratosCoSim.EMPIRE_API.EMPIRE_API_sendMesh(self.model[model_part_name], comm_name)

    def ImportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            variables = self.__GetVariables(interface_data.variable)
            KratosCoSim.EMPIRE_API.EMPIRE_API_recvDataField(interface_data.GetModelPart(), self.solver_name+"_"+interface_data.name, *variables)
        else:
            raise NotImplementedError('Importing interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def ExportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            variables = self.__GetVariables(interface_data.variable)
            KratosCoSim.EMPIRE_API.EMPIRE_API_sendDataField(interface_data.GetModelPart(), self.solver_name+"_"+interface_data.name, *variables)
        elif data_type == "convergence_signal":
            KratosCoSim.EMPIRE_API.EMPIRE_API_sendConvergenceSignal(data_config["is_converged"], self.solver_name)
        else:
            raise NotImplementedError('Exporting interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def __GetVariables(self, main_var):
        variables = [main_var]
        double_vector_var = self.double_vector_var_map.get(main_var)

        if double_vector_var is not None:
            variables.append(double_vector_var)

        return variables

    def PrintInfo(self):
        print("This is the EMPIRE-IO")

    def Check(self):
        pass

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "api_configuration_file_name" : "UNSPECIFIED",
            "communication_folder"        : ".EMPIRE"
            "communicate_as_double_vector : {}
        }""")
        this_defaults.AddMissingParameters(super(EmpireIO, cls)._GetDefaultSettings())

        return this_defaults
