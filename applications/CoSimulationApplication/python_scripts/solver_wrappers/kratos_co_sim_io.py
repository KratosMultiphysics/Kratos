# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_io import CoSimulationIO

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication import CoSimIO


def Create(settings, model, solver_name):
    return KratosCoSimIO(settings, model, solver_name)

class KratosCoSimIO(CoSimulationIO):
    """Wrapper for the CoSimIO to be used with Kratos
    """
    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)

        connect_to = self.settings["connect_to"].GetString()
        if connect_to == "":
            raise Exception('"connect_to" must be specified!')

        connection_settings = CoSimIO.InfoFromParameters(self.settings)
        connection_settings.SetString("my_name", solver_name)

        info = CoSimIO.Connect(connection_settings)
        if info.GetInt("connection_status") != CoSimIO.ConnectionStatus.Connected:
            raise Exception("Connecting failed!")

        self.connection_name = info.GetString("connection_name")

    def Finalize(self):
        disconnect_settings = CoSimIO.Info()
        disconnect_settings.SetString("connection_name", self.connection_name)

        info = CoSimIO.Disconnect(disconnect_settings)
        if info.GetInt("connection_status") != CoSimIO.ConnectionStatus.Disconnected:
            raise Exception("Disconnecting failed!")

    def ImportCouplingInterface(self, interface_config):
        model_part_name = interface_config["model_part_name"]

        info = CoSimIO.Info()
        info.SetString("connection_name", self.connection_name)
        info.SetString("identifier", model_part_name.replace(".", "-")) # TODO chec if better solution can be found

        CoSimIO.ImportMesh(info, self.model[model_part_name]) # TODO this can also be geometry at some point

    def ExportCouplingInterface(self, interface_config):
        model_part_name = interface_config["model_part_name"]

        info = CoSimIO.Info()
        info.SetString("connection_name", self.connection_name)
        info.SetString("identifier", model_part_name.replace(".", "-")) # TODO chec if better solution can be found

        CoSimIO.ExportMesh(info, self.model[model_part_name]) # TODO this can also be geometry at some point

    def ImportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            info = CoSimIO.Info()
            info.SetString("connection_name", self.connection_name)
            info.SetString("identifier", interface_data.name)

            CoSimIO.ImportData(info, interface_data.GetModelPart(), interface_data.variable, GetDataLocation(interface_data.location))
        else:
            raise NotImplementedError('Exporting interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def ExportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            info = CoSimIO.Info()
            info.SetString("connection_name", self.connection_name)
            info.SetString("identifier", interface_data.name)

            CoSimIO.ExportData(info, interface_data.GetModelPart(), interface_data.variable, GetDataLocation(interface_data.location))

        elif data_type == "control_signal":
            info = CoSimIO.Info()
            info.SetString("connection_name", self.connection_name)
            info.SetString("identifier", "run_control")
            info.SetString("control_signal", data_config["control_signal"])
            settings = data_config.get("settings")
            if settings:
                info.SetInfo("settings", CoSimIO.InfoFromParameters(settings))

            CoSimIO.ExportInfo(info)

        elif data_type == "repeat_time_step":
            info = CoSimIO.Info()
            info.SetString("connection_name", self.connection_name)
            info.SetString("identifier", "repeat_time_step_info")
            info.SetBool("repeat_time_step", data_config["repeat_time_step"])

            CoSimIO.ExportInfo(info)
        else:
            raise NotImplementedError('Exporting interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def PrintInfo(self):
        print("This is the KratosCoSimIO")

    def Check(self):
        pass

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "connect_to"           : "",
            "communication_format" : "file",
            "print_timing"         : false
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults

def GetDataLocation(location_str):
    location_map = {
        "node_historical"     : CoSimIO.DataLocation.NodeHistorical,
        "node_non_historical" : CoSimIO.DataLocation.NodeNonHistorical,
        "element"             : CoSimIO.DataLocation.Element,
        "condition"           : CoSimIO.DataLocation.Condition,
        "model_part"          : CoSimIO.DataLocation.ModelPart
    }
    return location_map[location_str]
