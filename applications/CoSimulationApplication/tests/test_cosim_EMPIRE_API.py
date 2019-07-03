from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.CoSimulationApplication as KratosCoSim

import os

class TestCoSim_EMPIRE_API(KratosUnittest.TestCase):

    def test_unused_fcts(self):
        KratosCoSim.EMPIRE_API.EMPIRE_API_Connect("dummy.xml")
        KratosCoSim.EMPIRE_API.EMPIRE_API_Disconnect()
        KratosCoSim.EMPIRE_API.EMPIRE_API_getUserDefinedText("dummy_element")

    def test_EMPIRE_API_sendConvergenceSignal(self):
        with self.assertRaisesRegex(RuntimeError, "Input can only be 0 non-convergence or 1 convergence, called with: 2"):
            KratosCoSim.EMPIRE_API.EMPIRE_API_sendConvergenceSignal(2)

        signal = 0 # means convergence
        conv_signal_file_name = "EMPIRE_convergence_signal.dat" # this is hardcoded in C++

        KratosCoSim.EMPIRE_API.EMPIRE_API_sendConvergenceSignal(signal)

        self.assertTrue(os.path.isfile(conv_signal_file_name))

        with open(conv_signal_file_name, 'r') as conv_signal_file:
            content = conv_signal_file.read()
            self.assertEqual(content, str(signal))

        os.remove(conv_signal_file_name)


if __name__ == '__main__':
    KratosUnittest.main()
