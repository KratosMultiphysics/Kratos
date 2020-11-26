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
        super(KratosCoSimIO, self).__init__(settings, model, solver_name)

        connection_settings = CoSimIO.InfoFromParameters(self.settings)
        connection_settings.SetString("connection_name", solver_name)

        info = CoSimIO.Connect(connection_settings)
        if info.GetInt("connection_status") != CoSimIO.ConnectionStatus.Connected:
            raise Exception("Connecting failed!")

    def Finalize(self):
        disconnect_settings = CoSimIO.Info()
        disconnect_settings.SetString("connection_name", self.solver_name)

        info = CoSimIO.Disconnect(disconnect_settings)
        if info.GetInt("connection_status") != CoSimIO.ConnectionStatus.Disconnected:
            raise Exception("Disconnecting failed!")

    def ImportCouplingInterface(self, interface_config):
        model_part_name = interface_config["model_part_name"]

        info = CoSimIO.Info()
        info.SetString("connection_name", self.solver_name)
        info.SetString("identifier", model_part_name)

        CoSimIO.ImportMesh(info, self.model[model_part_name]) # TODO this can also be geometry at some point

    def ExportCouplingInterface(self, interface_config):
        model_part_name = interface_config["model_part_name"]

        info = CoSimIO.Info()
        info.SetString("connection_name", self.solver_name)
        info.SetString("identifier", model_part_name)

        CoSimIO.ExportMesh(info, self.model[model_part_name]) # TODO this can also be geometry at some point

    def ImportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            info = CoSimIO.Info()
            info.SetString("connection_name", self.solver_name)
            info.SetString("identifier", interface_data.name)

            CoSimIO.ImportData(info, interface_data.GetModelPart(), interface_data.variable, GetDataLocation(interface_data.location))

        elif data_type == "time":
            time_list = CoSimIO.ImportData(self.solver_name, "time_to_co_sim")
            if len(time_list) != 1:
                raise Exception("Wrong size received!")
            data_config["time"] = time_list[0]
        else:
            raise NotImplementedError('Exporting interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def ExportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            info = CoSimIO.Info()
            info.SetString("connection_name", self.solver_name)
            info.SetString("identifier", interface_data.name)

            CoSimIO.ExportData(info, interface_data.GetModelPart(), interface_data.variable, GetDataLocation(interface_data.location))

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
            info = CoSimIO.Info()
            info.SetString("connection_name", self.solver_name)
            CoSimIO.SendControlSignal(info, control_signal_key)
        else:
            raise NotImplementedError('Exporting interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def PrintInfo(self):
        print("This is the KratosCoSimIO")

    def Check(self):
        pass

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "is_connection_master" : true,
            "communication_format" : "file",
            "print_timing"         : false
        }""")
        this_defaults.AddMissingParameters(super(KratosCoSimIO, cls)._GetDefaultParameters())
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
