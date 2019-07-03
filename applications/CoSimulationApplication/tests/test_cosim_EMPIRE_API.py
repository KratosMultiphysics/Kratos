from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.CoSimulationApplication as KratosCoSim

import os

conv_signal_file_name = "EMPIRE_convergence_signal.dat" # this is hardcoded in C++


class TestCoSim_EMPIRE_API(KratosUnittest.TestCase):

    def test_unused_fcts(self):
        KratosCoSim.EMPIRE_API.EMPIRE_API_Connect("dummy.xml")
        KratosCoSim.EMPIRE_API.EMPIRE_API_Disconnect()
        KratosCoSim.EMPIRE_API.EMPIRE_API_getUserDefinedText("dummy_element")

    def test_EMPIRE_API_sendConvergenceSignal(self):
        with self.assertRaisesRegex(RuntimeError, "Input can only be 0 for non-convergence or 1 for convergence, called with: 2"):
            KratosCoSim.EMPIRE_API.EMPIRE_API_sendConvergenceSignal(2)

        signal = 1 # means convergence

        KratosCoSim.EMPIRE_API.EMPIRE_API_sendConvergenceSignal(signal)

        self.__CheckConvergenceSignalFile(signal)

        os.remove(conv_signal_file_name)

    def test_EMPIRE_API_recvConvergenceSignal(self):
        # write file with convergence info manually
        signal = 0 # means non-convergence

        with open(conv_signal_file_name, 'w') as conv_signal_file:
            conv_signal_file.write(str(signal))

        self.__CheckConvergenceSignalFile(signal) # to make sure that the file was written correctly

        is_converged = KratosCoSim.EMPIRE_API.EMPIRE_API_recvConvergenceSignal()

        self.assertEqual(signal, is_converged)
        # make sure that the file was deleted
        self.assertFalse(os.path.isfile(conv_signal_file_name))

        # now check with an invalid signal
        invalid_signal = 15 # invalid

        with open(conv_signal_file_name, 'w') as conv_signal_file:
            conv_signal_file.write(str(invalid_signal))

        self.__CheckConvergenceSignalFile(invalid_signal) # to make sure that the file was written correctly

        with self.assertRaisesRegex(RuntimeError, "Read an invalid convergence signal: {}, can only be 0 for non-convergence or 1 for convergence".format(invalid_signal)):
            KratosCoSim.EMPIRE_API.EMPIRE_API_recvConvergenceSignal()

        # manually deleting since the call above throws
        os.remove(conv_signal_file_name)


    def __CheckConvergenceSignalFile(self, signal):
        self.assertTrue(os.path.isfile(conv_signal_file_name))

        with open(conv_signal_file_name, 'r') as conv_signal_file:
            content = conv_signal_file.read()
            self.assertEqual(content, str(signal))




if __name__ == '__main__':
    KratosUnittest.main()
