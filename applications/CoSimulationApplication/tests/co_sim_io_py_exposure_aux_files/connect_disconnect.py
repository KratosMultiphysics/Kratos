from KratosMultiphysics.CoSimulationApplication import CoSimIO

connection_settings = CoSimIO.Info()
connection_settings.SetString("my_name", "partner_b")
connection_settings.SetString("connect_to", "partner_a")
connection_settings.SetInt("echo_level", 2)
info = CoSimIO.Connect(connection_settings)
connection_name = info.GetString("connection_name")
if info.GetInt("connection_status") != CoSimIO.ConnectionStatus.Connected:
    raise Exception("Connecting failed")

disconnect_settings = CoSimIO.Info()
disconnect_settings.SetString("connection_name", connection_name)

info = CoSimIO.Disconnect(disconnect_settings)
if info.GetInt("connection_status") != CoSimIO.ConnectionStatus.Disconnected:
    raise Exception("Disconnecting failed")
