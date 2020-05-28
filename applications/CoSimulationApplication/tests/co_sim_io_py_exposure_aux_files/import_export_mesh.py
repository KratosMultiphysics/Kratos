import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication import CoSimIO

model = KM.Model()
model_part = model.CreateModelPart("mp_test")

connection_settings = CoSimIO.Info()
connection_settings.SetString("connection_name", "im_exp_mesh")
connection_settings.SetInt("echo_level", 0)
CoSimIO.Connect(connection_settings)

import_info = CoSimIO.Info()
import_info.SetString("connection_name", "im_exp_mesh")
import_info.SetString("identifier", "mesh_exchange_1")
CoSimIO.ImportMesh(import_info, model_part)

# print(model_part)

export_info = CoSimIO.Info()
export_info.SetString("connection_name", "im_exp_mesh")
export_info.SetString("identifier", "mesh_exchange_2")
CoSimIO.ExportMesh(export_info, model_part)

disconnect_settings = CoSimIO.Info()
disconnect_settings.SetString("connection_name", "im_exp_mesh")

CoSimIO.Disconnect(disconnect_settings)
