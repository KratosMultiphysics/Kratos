from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication import CoSimIO
from KratosMultiphysics import kratos_utilities as kratos_utils

import os
import threading

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestCoSimIOPyExposure(KratosUnittest.TestCase):

    def test_Connect_Disconnect(self):
        def Connect_Disconnect():
            connection_settings = CoSimIO.Info()
            connection_settings.SetString("connection_name", "abcc")
            connection_settings.SetInt("echo_level", 0)
            CoSimIO.Connect(connection_settings)

            disconnect_settings = CoSimIO.Info()
            disconnect_settings.SetString("connection_name", "abcc")

            CoSimIO.Disconnect(disconnect_settings)

        t = StartInThread(Connect_Disconnect)

        connection_settings = CoSimIO.Info()
        connection_settings.SetString("connection_name", "abcc")
        connection_settings.SetInt("echo_level", 0)
        CoSimIO.Connect(connection_settings)

        disconnect_settings = CoSimIO.Info()
        disconnect_settings.SetString("connection_name", "abcc")

        CoSimIO.Disconnect(disconnect_settings)

        t.join()

    def test_Export_Import_Data_raw_values(self):
        def Import_Export_Data():
            connection_settings = CoSimIO.Info()
            connection_settings.SetString("connection_name", "abcc")
            connection_settings.SetInt("echo_level", 0)
            CoSimIO.Connect(connection_settings)

            import_info = CoSimIO.Info()
            import_info.SetString("connection_name", "abcc")
            import_info.SetString("identifier", "raw_data")
            imported_values = CoSimIO.ImportData(import_info)

            export_info = CoSimIO.Info()
            export_info.SetString("connection_name", "abcc")
            export_info.SetString("identifier", "raw_data")
            CoSimIO.ExportData(export_info, imported_values)

            disconnect_settings = CoSimIO.Info()
            disconnect_settings.SetString("connection_name", "abcc")

            CoSimIO.Disconnect(disconnect_settings)

        # t = StartInThread(Import_Export_Data)

        connection_settings = CoSimIO.Info()
        connection_settings.SetString("connection_name", "abcc")
        connection_settings.SetInt("echo_level", 0)
        CoSimIO.Connect(connection_settings)

        values = [1.0, 2.5, 3.3, -9.4]

        export_info = CoSimIO.Info()
        export_info.SetString("connection_name", "abcc")
        export_info.SetString("identifier", "raw_data")
        CoSimIO.ExportData(export_info, values)

        import_info = CoSimIO.Info()
        import_info.SetString("connection_name", "abcc")
        import_info.SetString("identifier", "raw_data")
        imported_values = CoSimIO.ImportData(import_info)

        disconnect_settings = CoSimIO.Info()
        disconnect_settings.SetString("connection_name", "abcc")

        CoSimIO.Disconnect(disconnect_settings)

        t.join()


def StartInThread(target_fct):
    t = threading.Thread(
        target=target_fct
    )

    t.start()
    return t



if __name__ == '__main__':
    KratosUnittest.main()
