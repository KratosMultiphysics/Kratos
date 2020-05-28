import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication import CoSimIO
from KratosMultiphysics import kratos_utilities as kratos_utils

import os
import subprocess

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestCoSimIOPyExposure(KratosUnittest.TestCase):

    def test_Connect_Disconnect(self):
        connection_settings = CoSimIO.Info()
        connection_settings.SetString("connection_name", "c_d_test")
        connection_settings.SetInt("echo_level", 0)
        CoSimIO.Connect(connection_settings)

        RunPythonInSubProcess("connect_disconnect")

        disconnect_settings = CoSimIO.Info()
        disconnect_settings.SetString("connection_name", "c_d_test")

        CoSimIO.Disconnect(disconnect_settings)

    def test_Export_Import_Data_raw_values(self):
        connection_settings = CoSimIO.Info()
        connection_settings.SetString("connection_name", "im_exp_data")
        connection_settings.SetInt("echo_level", 0)
        CoSimIO.Connect(connection_settings)

        values = [1.0, 2.5, 3.3, -9.4]

        export_info = CoSimIO.Info()
        export_info.SetString("connection_name", "im_exp_data")
        export_info.SetString("identifier", "raw_data")
        CoSimIO.ExportData(export_info, values)

        RunPythonInSubProcess("import_export_data")

        import_info = CoSimIO.Info()
        import_info.SetString("connection_name", "im_exp_data")
        import_info.SetString("identifier", "return_raw_data")
        imported_values = CoSimIO.ImportData(import_info)

        disconnect_settings = CoSimIO.Info()
        disconnect_settings.SetString("connection_name", "im_exp_data")

        CoSimIO.Disconnect(disconnect_settings)

        # checking the values after disconnecting to avoid deadlock
        self.assertVectorAlmostEqual(KM.Vector(values), KM.Vector(imported_values))

def RunPythonInSubProcess(python_script_name):
    if not python_script_name.endswith(".py"):
        python_script_name += ".py"
    cmd_list = ["python3", os.path.join("co_sim_io_py_exposure_aux_files", python_script_name)]
    subprocess.run(cmd_list, check=True) # crashes the calling script too, otherwise the error is silent


if __name__ == '__main__':
    KratosUnittest.main()
