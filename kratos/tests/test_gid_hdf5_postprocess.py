import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics import gid_hdf5_postprocess

import h5py
import numpy as np
import os


class TestGiDHDF5PostProcess(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KratosMultiphysics.Model()
        self.model_part = self.model.CreateModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        coordinates = [
            (0.0, 0.0, 0.0),   # 1 corner
            (1.0, 0.0, 0.0),   # 2 corner
            (0.0, 1.0, 0.0),   # 3 corner
            (0.5, 0.0, 0.0),   # 4 mid 1-2
            (0.5, 0.5, 0.0),   # 5 mid 2-3
            (0.0, 0.5, 0.0),   # 6 mid 3-1
            (1.0, 1.0, 0.0),   # 7 corner
            (0.5, 1.0, 0.0),   # 8 mid 3-7
            (1.0, 0.5, 0.0),   # 9 mid 7-2
        ]
        self.coordinates = coordinates
        for node_id, (x, y, z) in enumerate(coordinates, start=1):
            self.model_part.CreateNewNode(node_id, x, y, z)
        property_ = self.model_part.GetProperties()[0]
        self.model_part.CreateNewElement("Element2D6N", 1, [1, 2, 3, 4, 5, 6], property_)
        self.model_part.CreateNewElement("Element2D6N", 2, [2, 3, 7, 5, 8, 9], property_)
        self.output_file = "gid_hdf5_postprocess_test_{}.h5".format(self._testMethodName)

    def tearDown(self):
        if os.path.exists(self.output_file):
            os.remove(self.output_file)

    def _CreateProcess(self):
        parameters = KratosMultiphysics.Parameters(
            """
        {
            "Parameters": {
                "model_part_name": "Main",
                "output_name": "__OUTPUT_FILE__",
                "postprocess_parameters": {
                    "result_file_configuration": {
                        "output_control_type": "step",
                        "output_interval": 1,
                        "nodal_results": ["PRESSURE", "DISPLACEMENT"],
                        "nodal_nonhistorical_results": ["DENSITY", "VELOCITY"],
                        "nodal_flags_results": ["BOUNDARY"]
                    }
                }
            }
        }
        """.replace("__OUTPUT_FILE__", self.output_file))
        process = gid_hdf5_postprocess.Factory(parameters, self.model)
        process.ExecuteInitialize()
        return process

    def _SetSolution(self, step):
        for node in self.model_part.Nodes:
            i = float(node.Id)
            node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 0, i * float(step))
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, KratosMultiphysics.Array3([i * float(step), 0.0, 0.0]))
            node.SetValue(KratosMultiphysics.DENSITY, 10.0 * i * float(step))
            node.SetValue(KratosMultiphysics.VELOCITY, KratosMultiphysics.Array3([2.0 * i * float(step), 0.0, 0.0]))
            if node.Id % 2 == 0:
                node.Set(KratosMultiphysics.BOUNDARY)

    def _RunTwoSteps(self, process):
        for step in (1, 2):
            self.model_part.ProcessInfo[KratosMultiphysics.STEP] = step
            self.model_part.ProcessInfo[KratosMultiphysics.TIME] = float(step)
            self.model_part.CloneTimeStep(float(step))
            self._SetSolution(step)
            process.PrintOutput()
        process.f.close()

    def _Attr(self, group, name):
        return group.attrs[name].decode()

    def _AssertMesh(self, f):
        mesh_1 = f["Meshes"]["1"]
        coords = mesh_1["Coordinates"]
        np.testing.assert_array_equal(np.asarray(coords["1"][:]), np.arange(1, 10))
        np.testing.assert_allclose(coords["2"][:], [c[0] for c in self.coordinates])
        np.testing.assert_allclose(coords["3"][:], [c[1] for c in self.coordinates])
        np.testing.assert_allclose(coords["4"][:], [c[2] for c in self.coordinates])
        np.testing.assert_array_equal(np.asarray(mesh_1["Elements"]["1"][:]), np.array([1, 2]))

    def test_historical_nodal_results(self):
        process = self._CreateProcess()
        self._RunTwoSteps(process)

        with h5py.File(self.output_file, "r") as f:
            self._AssertMesh(f)
            self.assertEqual(self._Attr(f, "GiD Post Results File"), "1.1")
            results = f["Results"]

            pressure_group_1 = results["1"]
            self.assertEqual(self._Attr(pressure_group_1, "Name"), "PRESSURE")
            self.assertEqual(self._Attr(pressure_group_1, "ResultType"), "Scalar")
            self.assertEqual(self._Attr(pressure_group_1, "NumComponents"), "1")
            self.assertEqual(self._Attr(pressure_group_1, "Step"), "1")
            np.testing.assert_array_equal(pressure_group_1["1"][:], np.arange(1, 10, dtype=process.int_type))
            np.testing.assert_allclose(pressure_group_1["2"][:], np.arange(1, 10))

            displacement_group_1 = results["2"]
            self.assertEqual(self._Attr(displacement_group_1, "Name"), "DISPLACEMENT")
            self.assertEqual(self._Attr(displacement_group_1, "ResultType"), "Vector")
            self.assertEqual(self._Attr(displacement_group_1, "NumComponents"), "3")
            self.assertEqual(self._Attr(displacement_group_1, "Component 1"), "DISPLACEMENT_X")
            np.testing.assert_allclose(displacement_group_1["2"][:], np.arange(1, 10))
            np.testing.assert_allclose(displacement_group_1["3"][:], np.zeros(9))

            pressure_group_2 = results["7"]
            self.assertEqual(self._Attr(pressure_group_2, "Name"), "PRESSURE")
            self.assertEqual(self._Attr(pressure_group_2, "Step"), "2")
            np.testing.assert_allclose(pressure_group_2["2"][:], np.arange(1, 10) * 2.0)

    def test_nonhistorical_nodal_results(self):
        process = self._CreateProcess()
        self._RunTwoSteps(process)

        with h5py.File(self.output_file, "r") as f:
            self._AssertMesh(f)
            results = f["Results"]

            density_group_1 = results["3"]
            self.assertEqual(self._Attr(density_group_1, "Name"), "DENSITY")
            self.assertEqual(self._Attr(density_group_1, "ResultType"), "Scalar")
            self.assertEqual(self._Attr(density_group_1, "Step"), "1")
            np.testing.assert_allclose(density_group_1["2"][:], np.arange(1, 10) * 10.0)

            velocity_group_1 = results["4"]
            self.assertEqual(self._Attr(velocity_group_1, "Name"), "VELOCITY")
            self.assertEqual(self._Attr(velocity_group_1, "ResultType"), "Vector")
            self.assertEqual(self._Attr(velocity_group_1, "Component 2"), "VELOCITY_Y")
            np.testing.assert_allclose(velocity_group_1["2"][:], np.arange(1, 10) * 2.0)
            np.testing.assert_allclose(velocity_group_1["3"][:], np.zeros(9))

            density_group_2 = results["9"]
            self.assertEqual(self._Attr(density_group_2, "Name"), "DENSITY")
            self.assertEqual(self._Attr(density_group_2, "Step"), "2")
            np.testing.assert_allclose(density_group_2["2"][:], np.arange(1, 10) * 20.0)

    def test_flags_nodal_results(self):
        process = self._CreateProcess()
        self._RunTwoSteps(process)

        with h5py.File(self.output_file, "r") as f:
            self._AssertMesh(f)
            results = f["Results"]

            flag_group_1 = results["5"]
            self.assertEqual(self._Attr(flag_group_1, "Name"), "BOUNDARY")
            self.assertEqual(self._Attr(flag_group_1, "ResultType"), "Scalar")
            self.assertEqual(self._Attr(flag_group_1, "NumComponents"), "1")
            self.assertEqual(self._Attr(flag_group_1, "Step"), "1")
            self.assertEqual(flag_group_1["2"].dtype.kind, "f")
            np.testing.assert_array_equal(flag_group_1["2"][:], np.array([-1, 1, -1, 1, -1, 1, -1, 1, -1]))

            flag_group_2 = results["11"]
            self.assertEqual(self._Attr(flag_group_2, "Name"), "BOUNDARY")
            self.assertEqual(self._Attr(flag_group_2, "Step"), "2")
            np.testing.assert_array_equal(flag_group_2["2"][:], np.array([-1, 1, -1, 1, -1, 1, -1, 1, -1]))

    def test_unknown_flag_raises(self):
        parameters = KratosMultiphysics.Parameters(
            """
        {
            "Parameters": {
                "model_part_name": "Main",
                "output_name": "__OUTPUT_FILE__",
                "postprocess_parameters": {
                    "result_file_configuration": {
                        "nodal_flags_results": ["NOT_A_FLAG"]
                    }
                }
            }
        }
        """.replace("__OUTPUT_FILE__", self.output_file))
        self.assertRaises(Exception, gid_hdf5_postprocess.Factory, parameters, self.model)


if __name__ == "__main__":
    KratosUnittest.main()
