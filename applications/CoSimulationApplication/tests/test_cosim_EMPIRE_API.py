from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.CoSimulationApplication as KratosCoSim

class TestCoSim_EMPIRE_API(KratosUnittest.TestCase):

    def test_connect_fcts(self):
        KratosCoSim.EMPIRE_API.EMPIRE_API_Connect("dummy.xml")
        KratosCoSim.EMPIRE_API.EMPIRE_API_Disconnect()

if __name__ == '__main__':
    KratosUnittest.main()
