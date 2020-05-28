from KratosMultiphysics.CoSimulationApplication import CoSimIO

connection_settings = CoSimIO.Info()
connection_settings.SetString("connection_name", "c_d_test")
connection_settings.SetInt("echo_level", 0)
CoSimIO.Connect(connection_settings)

disconnect_settings = CoSimIO.Info()
disconnect_settings.SetString("connection_name", "c_d_test")

CoSimIO.Disconnect(disconnect_settings)
