from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.CoSimulationApplication as KratosCoSim

import os, filecmp

conv_signal_file_name = "EMPIRE_convergence_signal.dat" # this is hardcoded in C++

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

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
        model = KM.Model()
        model_part = model.CreateModelPart("For_Sending")

        severity = KM.Logger.GetDefaultOutput().GetSeverity()
        KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING) # mute MP-IO

        model_part_io = KM.ModelPartIO(GetFilePath("generic_mdpa_files/Mok_CFD"))
        model_part_io.ReadModelPart(model_part)

        KM.Logger.GetDefaultOutput().SetSeverity(severity)

        KratosCoSim.EMPIRE_API.EMPIRE_API_sendMesh(model_part)

        # can directly use filecmp because there are no decimal-number issues
        self.assertTrue(filecmp.cmp(GetFilePath("reference_files/EMPIRE_mesh_For_Sending.vtk_ref"), "EMPIRE_mesh_For_Sending.vtk"))

    def test_EMPIRE_API_recvMesh(self):
        model = KM.Model()
        model_part = model.CreateModelPart("For_Receiving")
        model_part_ref = model.CreateModelPart("For_Checking")

        KratosCoSim.EMPIRE_API.EMPIRE_API_recvMesh(model_part)

        # reading reference ModelPart (with which the ref-vtk was created)
        severity = KM.Logger.GetDefaultOutput().GetSeverity()
        KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING) # mute MP-IO
        model_part_io = KM.ModelPartIO(GetFilePath("generic_mdpa_files/Mok_CFD"))
        model_part_io.ReadModelPart(model_part_ref)
        KM.Logger.GetDefaultOutput().SetSeverity(severity)

        self.assertEqual(model_part.NumberOfNodes(), model_part_ref.NumberOfNodes())
        self.assertEqual(model_part.NumberOfElements(), model_part_ref.NumberOfElements())

        self.__CompareNodes(model_part.Nodes, model_part_ref.Nodes())

        for elem, elem_ref in zip(model_part.Elements, model_part_ref.Elements()):
            self.assertEqual(elem.Id, elem_ref.Id)
            self.__CompareNodes(elem.GetNodes(), elem_ref.GetNodes())

    def test_EMPIRE_API_sendDataField(self):
        fct_to_test = KratosCoSim.EMPIRE_API.EMPIRE_API_sendDataField
        file_name_fct_ptr = GetDataFieldFileName

        array = [1.0, 2.5, 3.5, -99.11, -0.02, 555.5, -99.114, 0.02, 565.5, 10.0, 78.44]

        self.__TestArraySend(fct_to_test, file_name_fct_ptr, array)

    def test_EMPIRE_API_recvDataField(self):
        fct_to_test = KratosCoSim.EMPIRE_API.EMPIRE_API_recvDataField
        file_name_fct_ptr = GetDataFieldFileName

        array = [1.0, 2.5, 3.5, -99.11, -0.02, 555.5, -99.114, -99.11, -0.02, 555.5, 0.02, 565.5, 10.0, 78.44]

        self.__TestArrayReceive(fct_to_test, file_name_fct_ptr, array)

    def test_EMPIRE_API_sendSignal_double(self):
        fct_to_test = KratosCoSim.EMPIRE_API.EMPIRE_API_sendSignal_double
        file_name_fct_ptr = GetSignalFileName

        array = [1.0, 2.5, 3.5, -99.11, -0.02, 555.5, -99.114, 0.02, 565.5, 10.0, 78.44]

        self.__TestArraySend(fct_to_test, file_name_fct_ptr, array)

    def test_EMPIRE_API_recvSignal_double(self):
        fct_to_test = KratosCoSim.EMPIRE_API.EMPIRE_API_recvSignal_double
        file_name_fct_ptr = GetSignalFileName

        array = [1.0, 2.5, 3.5, -99.11, -0.02, 555.5, -99.114, -99.11, -0.02, 555.5, 0.02, 565.5, 10.0, 78.44]

        self.__TestArrayReceive(fct_to_test, file_name_fct_ptr, array)


    def __CheckConvergenceSignalFile(self, signal):
        self.assertTrue(os.path.isfile(conv_signal_file_name))

        with open(conv_signal_file_name, 'r') as conv_signal_file:
            content = conv_signal_file.read()
            self.assertEqual(content, str(signal))

    def __TestArraySend(self, fct_ptr_to_test, file_name_fct_ptr, array_to_test):
        array_name = "dummy_array_send"
        array_file_name = file_name_fct_ptr(array_name)

        fct_ptr_to_test(array_name, len(array_to_test), array_to_test)

        self.assertTrue(os.path.isfile(array_file_name))

        with open(array_file_name, 'r') as array_file:
            content = array_file.read()
            vals = [float(v) for v in content.split(' ')]
            for v, v_exp in zip(vals, array_to_test):
                self.assertAlmostEqual(v, v_exp)

        os.remove(array_file_name)

    def __TestArrayReceive(self, fct_ptr_to_test, file_name_fct_ptr, array_to_test):
        array_name = "dummy_array_recv"
        array_file_name = file_name_fct_ptr(array_name)

        with open(array_file_name, 'w') as array_file:
            for i_v, v in enumerate(array_to_test):
                array_file.write(str(v))
                # doing this extra to not have a trailing whitespace in the file
                if i_v < len(array_to_test)-1:
                    array_file.write(" ")

        array_to_receive = [0.0] * len(array_to_test)

        fct_ptr_to_test(array_name, len(array_to_receive), array_to_receive)

        # check the received signal
        for v, v_exp in zip(array_to_receive, array_to_test):
            self.assertAlmostEqual(v, v_exp)

        # make sure that the file was deleted
        self.assertFalse(os.path.isfile(array_file_name))

        with self.assertRaisesRegex(RuntimeError, "The size of the list has to be specified before, expected size of {}, current size: {}".format(len(array_to_test)+2, len(array_to_test))):
            fct_ptr_to_test(array_name, len(array_to_test)+2, array_to_test)

    def __CompareNodes(self, nodes, nodes_ref):
        for node, node_ref in zip(nodes, nodes_ref):
            self.assertEqual(node.Id, node_ref.Id)

            self.assertAlmostEqual(node.X(), node_ref.X())
            self.assertAlmostEqual(node.Y(), node_ref.Y())
            self.assertAlmostEqual(node.Z(), node_ref.Z())

            self.assertAlmostEqual(node.X0(), node_ref.X0())
            self.assertAlmostEqual(node.Y0(), node_ref.Y0())
            self.assertAlmostEqual(node.Z0(), node_ref.Z0())


def GetSignalFileName(signal_name):
    return "EMPIRE_signal_" + signal_name + ".dat" # this is hardcoded in C++

def GetDataFieldFileName(data_field_name):
    return "EMPIRE_datafield_" + data_field_name + ".dat" # this is hardcoded in C++


if __name__ == '__main__':
    KratosUnittest.main()
