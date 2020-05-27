from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication import CoSimIO
from KratosMultiphysics import kratos_utilities as kratos_utils

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestCoSimIOKratosPythonExposure(KratosUnittest.TestCase):

    def test_Connect_Disconnect(self):
        connection_settings = CoSimIO.Info()
        connection_settings.SetString("connection_name", "abcc")
        CoSimIO.Connect(connection_settings)

        disconnect_settings = CoSimIO.Info()
        disconnect_settings.SetString("connection_name", "abcc")

        CoSimIO.Disconnect(disconnect_settings)

if __name__ == '__main__':
    KratosUnittest.main()
