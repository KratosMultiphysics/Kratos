import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.testing.utilities import GetPython3Command

from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos_co_sim_io import Create as CreateKratosCoSimIO
from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData

import os
import subprocess


class TestKratosCoSimIO(KratosUnittest.TestCase):

    def test_Connect_Disconnect(self):
        p = self.__RunPythonInSubProcess("connect_disconnect")
        model = KM.Model()
        io_settings = KM.Parameters("""{
            "connect_to" : "partner_b"
        }""")
        solver_name = "partner_a" # aka "my_name" for the CoSimIO
        kratos_co_sim_io = CreateKratosCoSimIO(io_settings, model, solver_name, KM.Testing.GetDefaultDataCommunicator()) # this connects
        kratos_co_sim_io.Initialize()

        kratos_co_sim_io.Finalize() # this disconnects

        self.__CheckSubProcess(p)

    def test_Export_ImportCouplingData(self):
        p = self.__RunPythonInSubProcess("import_export_data")

        model = KM.Model()
        io_settings = KM.Parameters("""{
            "connect_to" : "impExp"
        }""")
        solver_name = "ExpImp" # aka "my_name" for the CoSimIO
        kratos_co_sim_io = CreateKratosCoSimIO(io_settings, model, solver_name, KM.Testing.GetDefaultDataCommunicator()) # this connects
        kratos_co_sim_io.Initialize()

        model_part = model.CreateModelPart("for_test")
        model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        model_part.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        for i in range(5):
            node = model_part.CreateNewNode(i+1, 0.0, 0.0, 0.0) # using same coord, doesn't matter for this test
            node.SetSolutionStepValue(KM.PRESSURE, i*1.7)

        settings_pres = KM.Parameters("""{
            "model_part_name" : "for_test",
            "variable_name"   : "PRESSURE"
        }""")

        settings_temp = KM.Parameters("""{
            "model_part_name" : "for_test",
            "variable_name"   : "TEMPERATURE"
        }""")
        interface_data_pres = CouplingInterfaceData(settings_pres, model, "data_exchange_1")
        interface_data_temp = CouplingInterfaceData(settings_temp, model, "data_exchange_2")

        data_configuration_export = {"type" : "coupling_interface_data", "interface_data" : interface_data_pres}
        kratos_co_sim_io.ExportData(data_configuration_export)

        data_configuration_import = {"type" : "coupling_interface_data", "interface_data" : interface_data_temp}
        kratos_co_sim_io.ImportData(data_configuration_import)

        kratos_co_sim_io.Finalize() # this disconnects

        # checking the values after disconnecting to avoid deadlock
        for i, node in enumerate(model_part.Nodes):
            self.assertAlmostEqual(node.GetSolutionStepValue(KM.TEMPERATURE), i*1.7)

        self.__CheckSubProcess(p)

    def test_Export_ImportCouplingInterface_Mesh(self): # can also be Geometry at some point
        p = self.__RunPythonInSubProcess("import_export_mesh")

        model = KM.Model()
        io_settings = KM.Parameters("""{
            "connect_to" : "impExpMesh"
        }""")
        solver_name = "ExpImpMesh" # aka "my_name" for the CoSimIO
        kratos_co_sim_io = CreateKratosCoSimIO(io_settings, model, solver_name, KM.Testing.GetDefaultDataCommunicator()) # this connects
        kratos_co_sim_io.Initialize()

        model_part = model.CreateModelPart("mesh_exchange_1")
        model_part_returned = model.CreateModelPart("mesh_exchange_2")

        for i in range(15):
            model_part.CreateNewNode(i+1, i*1.1, 1.1**i, 0.0)
        props = model_part.CreateNewProperties(0)
        for i in range(9): # this leaves some hanging nodes
            model_part.CreateNewElement("Element2D2N", i+1, [i+1, i+2], props)

        interface_configuration_export = {"model_part_name" : "mesh_exchange_1"}
        kratos_co_sim_io.ExportCouplingInterface(interface_configuration_export)

        interface_configuration_import = {"model_part_name" : "mesh_exchange_2"}
        kratos_co_sim_io.ImportCouplingInterface(interface_configuration_import)

        kratos_co_sim_io.Finalize() # this disconnects

        # checking the values after disconnecting to avoid deadlock
        self.assertEqual(model_part.NumberOfNodes(), model_part_returned.NumberOfNodes())
        self.assertEqual(model_part.NumberOfElements(), model_part_returned.NumberOfElements())

        for node_orig, node_exchanged in zip(model_part.Nodes, model_part_returned.Nodes):
            self.assertEqual(node_orig.Id, node_exchanged.Id)
            self.assertAlmostEqual(node_orig.X0, node_exchanged.X0)
            self.assertAlmostEqual(node_orig.Y0, node_exchanged.Y0)
            self.assertAlmostEqual(node_orig.Z0, node_exchanged.Z0)

        for elem_orig, elem_exchanged in zip(model_part.Elements, model_part_returned.Elements):
            self.assertEqual(len(elem_orig.GetNodes()), len(elem_exchanged.GetNodes()))
            for node_orig, node_exchanged in zip(elem_orig.GetNodes(), elem_exchanged.GetNodes()):
                self.assertEqual(node_orig.Id, node_exchanged.Id)
                self.assertAlmostEqual(node_orig.X0, node_exchanged.X0)
                self.assertAlmostEqual(node_orig.Y0, node_exchanged.Y0)
                self.assertAlmostEqual(node_orig.Z0, node_exchanged.Z0)

        self.__CheckSubProcess(p)


    def __RunPythonInSubProcess(self, script_name):
        if not script_name.endswith(".py"):
            script_name += ".py"

        return subprocess.Popen([GetPython3Command(), os.path.join("co_sim_io_py_exposure_aux_files", script_name)], stdout=subprocess.PIPE)

    def __CheckSubProcess(self, proc):
        try:
            p_out = proc.communicate(timeout=5)
        except subprocess.TimeoutExpired: # Timeout reached
            proc.kill()
            p_out = proc.communicate()

        self.assertEqual(proc.returncode, 0, msg=p_out[0].decode('ascii'))


if __name__ == '__main__':
    KratosUnittest.main()
