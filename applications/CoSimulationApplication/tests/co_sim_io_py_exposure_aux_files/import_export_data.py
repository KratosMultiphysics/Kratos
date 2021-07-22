from KratosMultiphysics.CoSimulationApplication import CoSimIO

connection_settings = CoSimIO.Info()
connection_settings.SetString("my_name", "impExp")
connection_settings.SetString("connect_to", "ExpImp")
connection_settings.SetInt("echo_level", 0)
info = CoSimIO.Connect(connection_settings)
connection_name = info.GetString("connection_name")
if info.GetInt("connection_status") != CoSimIO.ConnectionStatus.Connected:
    raise Exception("Connecting failed")

import_info = CoSimIO.Info()
import_info.SetString("connection_name", connection_name)
import_info.SetString("identifier", "data_exchange_1")
imported_values = CoSimIO.DoubleVector()
CoSimIO.ImportData(import_info, imported_values)

# print(imported_values)

export_info = CoSimIO.Info()
export_info.SetString("connection_name", connection_name)
export_info.SetString("identifier", "data_exchange_2")
CoSimIO.ExportData(export_info, imported_values)

disconnect_settings = CoSimIO.Info()
disconnect_settings.SetString("connection_name", connection_name)

info = CoSimIO.Disconnect(disconnect_settings)
if info.GetInt("connection_status") != CoSimIO.ConnectionStatus.Disconnected:
    raise Exception("Disconnecting failed")
