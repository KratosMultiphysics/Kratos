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

    def test_EMPIRE_API_sendMesh(self):
        pass

    def test_EMPIRE_API_recvMesh(self):
        pass

    def test_EMPIRE_API_sendDataField(self):
        pass

    def test_EMPIRE_API_recvDataField(self):
        pass

    def test_EMPIRE_API_sendSignal_double(self):
        signal_name = "dummy_signal_send"
        signal_file_name = GetSignalFileName(signal_name)

        signal = [1.0, 2.5, 3.5, -99.11, -0.02, 555.5]
        KratosCoSim.EMPIRE_API.EMPIRE_API_sendSignal_double(signal_name, len(signal), signal)


        self.assertTrue(os.path.isfile(signal_file_name))

        with open(signal_file_name, 'r') as signal_file:
            content = signal_file.read()
            vals = [float(v) for v in content.split(' ')]
            for v, v_exp in zip(vals, signal):
                self.assertAlmostEqual(v, v_exp)

        os.remove(signal_file_name)


    def test_EMPIRE_API_recvSignal_double(self):
        signal_name = "dummy_signal_recv"
        signal_file_name = GetSignalFileName(signal_name)

        exp_signal = [13.0, -21.5, 3.555, -99.114, 0.02, 565.5, 10.0, 78.44]

        with open(signal_file_name, 'w') as signal_file:
            for i_v, v in enumerate(exp_signal):
                signal_file.write(str(v))
                # doing this extra to not have a trailing whitespace in the file
                if i_v < len(exp_signal)-1:
                    signal_file.write(" ")

        signal = [0.0] * len(exp_signal)

        KratosCoSim.EMPIRE_API.EMPIRE_API_recvSignal_double(signal_name, len(signal), signal)

        # check the received signal
        for v, v_exp in zip(signal, exp_signal):
            self.assertAlmostEqual(v, v_exp)

        # make sure that the file was deleted
        self.assertFalse(os.path.isfile(signal_file_name))

    def __CheckConvergenceSignalFile(self, signal):
        self.assertTrue(os.path.isfile(conv_signal_file_name))

        with open(conv_signal_file_name, 'r') as conv_signal_file:
            content = conv_signal_file.read()
            self.assertEqual(content, str(signal))

def GetSignalFileName(signal_name):
    return "EMPIRE_signal_" + signal_name # this is hardcoded in C++

if __name__ == '__main__':
    KratosUnittest.main()
