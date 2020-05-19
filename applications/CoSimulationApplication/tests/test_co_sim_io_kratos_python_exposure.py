from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

import KratosMultiphysics.CoSimulationApplication as KratosCoSim
from KratosMultiphysics import kratos_utilities as kratos_utils

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestCoSimIOKratosPythonExposure(KratosUnittest.TestCase):

    # TODO add more tests, similar to the EMPIRE_API-tests

    def test_Connect_Disconnect(self):
        connection_name = "abcc"
        KratosCoSim.CoSimIO.Connect(connection_name, {"echo_level" : "0"})
        KratosCoSim.CoSimIO.Disconnect(connection_name)

if __name__ == '__main__':
    KratosUnittest.main()
