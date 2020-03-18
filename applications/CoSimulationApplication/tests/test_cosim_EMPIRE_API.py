from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

import KratosMultiphysics.CoSimulationApplication as KratosCoSim
from KratosMultiphysics import kratos_utilities as kratos_utils

import os
from shutil import copyfile

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

communication_folder = ".EmpireIO" # hardcoded in C++
conv_signal_file_name = os.path.join(communication_folder, "EMPIRE_convergence_signal_default.dat")

class TestCoSim_EMPIRE_API(KratosUnittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # to silence prints
        KratosCoSim.EMPIRE_API.EMPIRE_API_SetEchoLevel(0)
        KratosCoSim.EMPIRE_API.EMPIRE_API_PrintTiming(False)

    def setUp(self):
        # delete and recreate communication folder to avoid leftover files
        kratos_utils.DeleteDirectoryIfExisting(communication_folder)
        os.mkdir(communication_folder)

    def tearDown(self):
        kratos_utils.DeleteDirectoryIfExisting(communication_folder)


    def test_unused_fcts(self):
        # to make sure theses fcts still exist in case a solver still calls them
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

        params = KM.Parameters("""{
            "comparison_type"     : "vtk",
            "reference_file_name" : "",
            "output_file_name"    : ""
        }""")
        params["reference_file_name"].SetString(GetFilePath("reference_files/EMPIRE_mesh_For_Sending.vtk_ref"))
        params["output_file_name"].SetString(os.path.join(communication_folder, "EMPIRE_mesh_For_Sending.vtk"))

        CompareTwoFilesCheckProcess(params).Execute()

        kratos_utils.DeleteFileIfExisting("EMPIRE_mesh_For_Sending.vtk")

    def test_EMPIRE_API_recvMesh_std_vector(self):
        self.__TestReceiveMesh(False)

    def test_EMPIRE_API_recvMesh_raw_pointers(self):
        self.__TestReceiveMesh(True)

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

    def test_SendRecv_scalar_datafield(self):
        model = KM.Model()
        model_part_send = model.CreateModelPart("For_Sending")
        model_part_recv = model.CreateModelPart("For_Receiving")
        FillModelPart(model_part_send)
        InitializeModelPart(model_part_recv)

        KratosCoSim.EMPIRE_API.EMPIRE_API_sendDataField(model_part_send, KM.PRESSURE)
        KratosCoSim.EMPIRE_API.EMPIRE_API_recvDataField(model_part_recv, KM.PRESSURE)

        self.__CompareScalarNodalValues(model_part_recv.Nodes, model_part_send.Nodes, KM.PRESSURE)

    def test_SendRecv_vector_datafield(self):
        model = KM.Model()
        model_part_send = model.CreateModelPart("For_Sending")
        model_part_recv = model.CreateModelPart("For_Receiving")
        FillModelPart(model_part_send)
        InitializeModelPart(model_part_recv)

        KratosCoSim.EMPIRE_API.EMPIRE_API_sendDataField(model_part_send, KM.DISPLACEMENT)
        KratosCoSim.EMPIRE_API.EMPIRE_API_recvDataField(model_part_recv, KM.DISPLACEMENT)

        self.__CompareVectorNodalValues(model_part_recv.Nodes, model_part_send.Nodes, KM.DISPLACEMENT)

    def test_SendRecv_doubleVector_datafield(self):
        model = KM.Model()
        model_part_send = model.CreateModelPart("For_Sending")
        model_part_recv = model.CreateModelPart("For_Receiving")
        FillModelPart(model_part_send)
        InitializeModelPart(model_part_recv)

        KratosCoSim.EMPIRE_API.EMPIRE_API_sendDataField(model_part_send, KM.DISPLACEMENT, KM.ROTATION)
        KratosCoSim.EMPIRE_API.EMPIRE_API_recvDataField(model_part_recv, KM.DISPLACEMENT, KM.ROTATION)

        self.__CompareVectorNodalValues(model_part_recv.Nodes, model_part_send.Nodes, KM.DISPLACEMENT)
        self.__CompareVectorNodalValues(model_part_recv.Nodes, model_part_send.Nodes, KM.ROTATION)


    def __CheckConvergenceSignalFile(self, signal):
        self.assertTrue(os.path.isfile(conv_signal_file_name))

        with open(conv_signal_file_name, 'r') as conv_signal_file:
            content = conv_signal_file.read()
            self.assertEqual(content, str(signal))

    def __TestArraySend(self, fct_ptr_to_test, file_name_fct_ptr, array_to_test):
        array_name = "dummy_array_send"
        array_file_name = file_name_fct_ptr(array_name)
        array_size = len(array_to_test)

        fct_ptr_to_test(array_name, array_size, array_to_test)

        self.assertTrue(os.path.isfile(array_file_name))

        with open(array_file_name, 'r') as array_file:
            content = array_file.readlines()
            self.assertEqual(len(content), 2) # file should contain 2 lines
            written_array_size = int(content[0])
            self.assertEqual(array_size, written_array_size)
            vals = [float(v) for v in content[1].split(' ')]
            self.assertEqual(array_size, len(vals))
            for v, v_exp in zip(vals, array_to_test):
                self.assertAlmostEqual(v, v_exp, 10)

        os.remove(array_file_name)

    def __TestArrayReceive(self, fct_ptr_to_test, file_name_fct_ptr, array_to_test):
        array_name = "dummy_array_recv"
        array_file_name = file_name_fct_ptr(array_name)

        array_size = len(array_to_test)
        with open(array_file_name, 'w') as array_file:
            array_file.write(str(array_size)+"\n")
            for i_v, v in enumerate(array_to_test):
                array_file.write(str(v))
                # doing this extra to not have a trailing whitespace in the file
                if i_v < array_size-1:
                    array_file.write(" ")

        array_to_receive = [0.0] * array_size
        fct_ptr_to_test(array_name, len(array_to_receive), array_to_receive)

        # check the received signal
        for v, v_exp in zip(array_to_receive, array_to_test):
            self.assertAlmostEqual(v, v_exp, 10)

        # make sure that the file was deleted
        self.assertFalse(os.path.isfile(array_file_name))

        # writting a file with the wrong size, this should throw a proper error
        with open(array_file_name, 'w') as array_file:
            array_file.write(str(array_size+2)+"\n")
        with self.assertRaisesRegex(RuntimeError, 'The received size for array "{}" is different from what is expected:\n    Expected size: {}\n    Received size: {}'.format(array_file_name, array_size, array_size+2)):
            fct_ptr_to_test(array_name, len(array_to_receive), array_to_receive)

        # manually remove file since the call above throws
        os.remove(array_file_name)

        # check if the py-exposure works correctly
        with self.assertRaisesRegex(RuntimeError, "The size of the list has to be specified before, expected size of {}, current size: {}".format(array_size+2, array_size)):
            fct_ptr_to_test(array_name, array_size+2, array_to_test)

    def __TestReceiveMesh(self, use_raw_pointers):
        mp_name = "For_Receiving"
        mesh_file_name = GetMeshFileName(mp_name)
        model = KM.Model()
        model_part = model.CreateModelPart(mp_name)
        model_part_ref = model.CreateModelPart("For_Checking")

        copyfile(GetFilePath("reference_files/EMPIRE_mesh_For_Sending.vtk_ref"), mesh_file_name)

        KratosCoSim.EMPIRE_API.EMPIRE_API_recvMesh(model_part, False, use_raw_pointers)

        # make sure that the file was deleted
        self.assertFalse(os.path.isfile(mesh_file_name))

        # reading reference ModelPart (with which the ref-vtk was created)
        severity = KM.Logger.GetDefaultOutput().GetSeverity()
        KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING) # mute MP-IO
        model_part_io = KM.ModelPartIO(GetFilePath("generic_mdpa_files/Mok_CFD"))
        model_part_io.ReadModelPart(model_part_ref)
        KM.Logger.GetDefaultOutput().SetSeverity(severity)

        self.assertEqual(model_part.NumberOfNodes(), model_part_ref.NumberOfNodes())
        self.assertEqual(model_part.NumberOfElements(), model_part_ref.NumberOfElements())

        self.__CompareNodes(model_part.Nodes, model_part_ref.Nodes)

        for elem, elem_ref in zip(model_part.Elements, model_part_ref.Elements):
            self.assertEqual(elem.Id, elem_ref.Id)
            self.__CompareNodes(elem.GetNodes(), elem_ref.GetNodes())

    def __CompareNodes(self, nodes, nodes_ref):
        for node, node_ref in zip(nodes, nodes_ref):
            self.assertEqual(node.Id, node_ref.Id)

            self.assertAlmostEqual(node.X, node_ref.X)
            self.assertAlmostEqual(node.Y, node_ref.Y)
            self.assertAlmostEqual(node.Z, node_ref.Z)

            self.assertAlmostEqual(node.X0, node_ref.X0)
            self.assertAlmostEqual(node.Y0, node_ref.Y0)
            self.assertAlmostEqual(node.Z0, node_ref.Z0)

    def __CompareScalarNodalValues(self, nodes, nodes_ref, var):
        for node, node_ref in zip(nodes, nodes_ref):
            self.assertAlmostEqual(node.GetSolutionStepValue(var), node_ref.GetSolutionStepValue(var), 12)

    def __CompareVectorNodalValues(self, nodes, nodes_ref, var):
        for node, node_ref in zip(nodes, nodes_ref):
            val = node.GetSolutionStepValue(var)
            val_ref = node_ref.GetSolutionStepValue(var)
            for v, v_ref in zip(val, val_ref):
                self.assertAlmostEqual(v, v_ref, 12)


def GetPRESUREValue(node_id):
    return node_id * 10.458

def GetDISPLACEMENTValue(node_id):
    return [node_id*0.0001458, node_id+6, node_id-8569]

def GetROTATIONValue(node_id):
    return [node_id*0.00000005561458, node_id+613.9, node_id-0.0008569]

def InitializeModelPart(model_part):
    model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
    model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KM.ROTATION)

    for i in range(50):
        model_part.CreateNewNode(i+1, i*1.1, i*1.2, i*1.3)

def FillModelPart(model_part):
    InitializeModelPart(model_part)

    for node in model_part.Nodes:
        node_id = node.Id
        node.SetSolutionStepValue(KM.PRESSURE, GetPRESUREValue(node_id))
        node.SetSolutionStepValue(KM.DISPLACEMENT, GetDISPLACEMENTValue(node_id))
        node.SetSolutionStepValue(KM.ROTATION, GetROTATIONValue(node_id))

def GetSignalFileName(signal_name):
    return os.path.join(communication_folder, "EMPIRE_signal_" + signal_name + ".dat") # this is hardcoded in C++

def GetDataFieldFileName(data_field_name):
    return os.path.join(communication_folder, "EMPIRE_datafield_" + data_field_name + ".dat") # this is hardcoded in C++

def GetMeshFileName(mesh_name):
    return os.path.join(communication_folder, "EMPIRE_mesh_" + mesh_name + ".vtk") # this is hardcoded in C++


if __name__ == '__main__':
    KratosUnittest.main()
