import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos_co_sim_io import Create as CreateKratosCoSimIO
from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData

import os
import subprocess
from shutil import which


class TestKratosCoSimIO(KratosUnittest.TestCase):

    def test_Connect_Disconnect(self):
        model = KM.Model()
        io_settings = KM.Parameters("""{}""") # using the defaults
        solver_name = "c_d_test"
        kratos_co_sim_io = CreateKratosCoSimIO(io_settings, model, solver_name) # this connects
        kratos_co_sim_io.Initialize()

        RunPythonInSubProcess("connect_disconnect")

        kratos_co_sim_io.Finalize() # this disconnects

    def test_Export_ImportCouplingData(self):
        model = KM.Model()
        io_settings = KM.Parameters("""{}""") # using the defaults
        solver_name = "im_exp_data"
        kratos_co_sim_io = CreateKratosCoSimIO(io_settings, model, solver_name) # this connects
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
        interface_data_pres.Initialize()
        interface_data_temp.Initialize()

        data_configuration_export = {"type" : "coupling_interface_data", "interface_data" : interface_data_pres}
        kratos_co_sim_io.ExportData(data_configuration_export)

        RunPythonInSubProcess("import_export_data")

        data_configuration_import = {"type" : "coupling_interface_data", "interface_data" : interface_data_temp}
        kratos_co_sim_io.ImportData(data_configuration_import)

        kratos_co_sim_io.Finalize() # this disconnects

        # checking the values after disconnecting to avoid deadlock
        for i, node in enumerate(model_part.Nodes):
            self.assertAlmostEqual(node.GetSolutionStepValue(KM.TEMPERATURE), i*1.7)

    def test_Export_ImportCouplingInterface_Mesh(self): # can also be Geometry at some point
        model = KM.Model()
        io_settings = KM.Parameters("""{}""") # using the defaults
        solver_name = "im_exp_mesh"
        kratos_co_sim_io = CreateKratosCoSimIO(io_settings, model, solver_name) # this connects
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

        RunPythonInSubProcess("import_export_mesh")

        interface_configuration_import = {"model_part_name" : "mesh_exchange_2"}
        kratos_co_sim_io.ImportCouplingInterface(interface_configuration_import)

        kratos_co_sim_io.Finalize() # this disconnects

        # checking the values after disconnecting to avoid deadlock
        self.assertEqual(model_part.NumberOfNodes(), model_part_returned.NumberOfNodes())
        self.assertEqual(model_part.NumberOfElements(), model_part_returned.NumberOfElements())

        for node_orig, node_exchanged in zip(model_part.Nodes, model_part_returned.Nodes):
            self.assertAlmostEqual(node_orig.X0, node_exchanged.X0)
            self.assertAlmostEqual(node_orig.Y0, node_exchanged.Y0)
            self.assertAlmostEqual(node_orig.Z0, node_exchanged.Z0)

        for elem_orig, elem_exchanged in zip(model_part.Elements, model_part_returned.Elements):
            self.assertEqual(len(elem_orig.GetNodes()), len(elem_exchanged.GetNodes()))
            for node_orig, node_exchanged in zip(elem_orig.GetNodes(), elem_exchanged.GetNodes()):
                self.assertAlmostEqual(node_orig.X0, node_exchanged.X0)
                self.assertAlmostEqual(node_orig.Y0, node_exchanged.Y0)
                self.assertAlmostEqual(node_orig.Z0, node_exchanged.Z0)


def RunPythonInSubProcess(python_script_name):
    if not python_script_name.endswith(".py"):
        python_script_name += ".py"

    py_cmd = "python3" if which("python3") is not None else "python"

    cmd_list = [py_cmd, os.path.join("co_sim_io_py_exposure_aux_files", python_script_name)]
    subprocess.run(cmd_list, check=True, shell=os.name=="nt") # crashes the calling script too, otherwise the error is silent (using shell in Win)


if __name__ == '__main__':
    KratosUnittest.main()
