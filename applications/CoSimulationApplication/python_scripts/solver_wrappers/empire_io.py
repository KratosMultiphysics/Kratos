# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics import kratos_utilities

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_io import CoSimulationIO

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication as KratosCoSim
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities

# Other imports
import os

def Create(*args):
    return EmpireIO(*args)

communication_folder = ".EmpireIO" # hardcoded in C++

class EmpireIO(CoSimulationIO):
    """IO for the legacy EMPIRE_API
    """
    def __init__(self, settings, model, solver_name, data_commnicator):
        super().__init__(settings, model, solver_name, data_commnicator)
        # Note: calling "EMPIRE_API_Connect" is NOT necessary, it is replaced by the next two lines
        KratosCoSim.EMPIRE_API.EMPIRE_API_SetEchoLevel(self.echo_level)
        KratosCoSim.EMPIRE_API.EMPIRE_API_PrintTiming(self.settings["api_print_timing"].GetBool())

        # delete and recreate communication folder to avoid leftover files
        kratos_utilities.DeleteDirectoryIfExisting(communication_folder)
        os.mkdir(communication_folder)

    def Finalize(self):
        kratos_utilities.DeleteDirectoryIfExisting(communication_folder)

    def __del__(self):
        # make sure no communication files are left even if simulation is terminated prematurely
        if os.path.isdir(communication_folder):
            kratos_utilities.DeleteDirectoryIfExisting(communication_folder)
            if self.echo_level > 0:
                cs_tools.cs_print_info(self._ClassName(), "Deleting Communication folder in destructor")

    def ImportCouplingInterface(self, interface_config):
        model_part_name = interface_config["model_part_name"]
        comm_name = interface_config["comm_name"]

        if not self.model.HasModelPart(model_part_name):
            main_model_part_name, *sub_model_part_names = model_part_name.split(".")
            model_part_utilities.RecursiveCreateModelParts(self.model[main_model_part_name], ".".join(sub_model_part_names))

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
            KratosCoSim.EMPIRE_API.EMPIRE_API_recvDataField(interface_data.GetModelPart(), self.solver_name+"_"+interface_data.name, interface_data.variable)
        else:
            raise NotImplementedError('Importing interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def ExportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            KratosCoSim.EMPIRE_API.EMPIRE_API_sendDataField(interface_data.GetModelPart(), self.solver_name+"_"+interface_data.name, interface_data.variable)
        elif data_type == "repeat_time_step":
            KratosCoSim.EMPIRE_API.EMPIRE_API_sendConvergenceSignal((not data_config["repeat_time_step"]), self.solver_name)
        else:
            raise NotImplementedError('Exporting interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def PrintInfo(self):
        print("This is the EMPIRE-IO")

    def Check(self):
        pass

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "api_print_timing" : false
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())

        return this_defaults
