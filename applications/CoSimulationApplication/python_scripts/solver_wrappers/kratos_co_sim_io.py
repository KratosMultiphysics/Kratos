# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_io import CoSimulationIO

# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication import CoSimIO


def Create(*args):
    return KratosCoSimIO(*args)

class KratosCoSimIO(CoSimulationIO):
    """Wrapper for the CoSimIO to be used with Kratos
    """
    def __init__(self, settings, model, solver_name, data_communicator):
        # backward compatibility
        for param in ("connect_to", "communication_format", "print_timing"):
            if settings.Has(param):
                if not settings.Has("co_sim_io_settings"):
                    settings.AddEmptyValue("co_sim_io_settings")
                co_sim_io_settings = settings["co_sim_io_settings"]
                co_sim_io_settings.AddValue(param, settings[param])
                settings.RemoveValue(param)

        super().__init__(settings, model, solver_name, data_communicator)

        co_sim_io_settings = settings["co_sim_io_settings"]

        if not co_sim_io_settings.Has("my_name"):
            co_sim_io_settings.AddEmptyValue("my_name").SetString(solver_name)

        connection_settings = CoSimIO.InfoFromParameters(co_sim_io_settings)

        if self.data_communicator.IsDistributed():
            from KratosMultiphysics.CoSimulationApplication.MPIExtension import CoSimIO as CoSimIOMPI
            info = CoSimIOMPI.ConnectMPI(connection_settings, self.data_communicator)
        else:
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

        CoSimIO.ImportMesh(info, self.model[model_part_name], self.data_communicator) # TODO this can also be geometry at some point

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
            "co_sim_io_settings" : { }
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
