import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.MedApplication as KratosMed

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from testing_utilities import MedModelPartIOTestCase, GetMedPath, get_num_geometries_by_type

from pathlib import Path


class TestMedModelPartIO(MedModelPartIOTestCase):
    def setUp(self):
        self.model = KM.Model()
        self.mp_read_1 = self.model.CreateModelPart("read_1")
        self.mp_read_2 = self.model.CreateModelPart("read_2")

    def _execute_tests(self, med_path, check_fct, print_vtk=False):
        med_io_read_1 = KratosMed.MedModelPartIO(GetMedPath(med_path))
        med_io_read_1.ReadModelPart(self.mp_read_1)

        self._basic_checks(self.mp_read_1)

        if print_vtk:
            write_vtk(self.mp_read_1, med_path)

        with self.subTest("check_model_part"):
            check_fct(self.mp_read_1)

        with self.subTest("read_write_read"):
            med_temp_path = GetMedPath(med_path).with_suffix(".med.tmp")
            DeleteFileIfExisting(med_temp_path)  # make sure there are no leftovers from previous tests
            self.addCleanup(DeleteFileIfExisting, med_temp_path)  # clean up after test

            med_io_write = KratosMed.MedModelPartIO(med_temp_path, KM.IO.WRITE)
            med_io_write.WriteModelPart(self.mp_read_1)

            med_io_read_2 = KratosMed.MedModelPartIO(med_temp_path, KM.IO.READ)
            med_io_read_2.ReadModelPart(self.mp_read_2)

            KratosMed.MedTestingUtilities.CheckModelPartsAreEqual(
                self.mp_read_1, self.mp_read_2, check_sub_model_parts=False  # until writing of SMPs is implemented
            )

    def test_empty_med_file(self):
        def mp_check_fct(model_part):
            self.assertEqual(model_part.NumberOfNodes(), 0)

        self._execute_tests(Path("empty"), mp_check_fct)

    def test_only_nodes(self):
        def mp_check_fct(model_part):
            self.assertEqual(model_part.NumberOfNodes(), 4)

            exp_coords = [(0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 1, 1)]

            for coords, node in zip(exp_coords, model_part.Nodes):
                self.assertAlmostEqual(node.X, coords[0])
                self.assertAlmostEqual(node.X0, coords[0])
                self.assertAlmostEqual(node.Y, coords[1])
                self.assertAlmostEqual(node.Y0, coords[1])
                self.assertAlmostEqual(node.Z, coords[2])
                self.assertAlmostEqual(node.Z0, coords[2])

        self._execute_tests(Path("only_nodes"), mp_check_fct)

    def test_nodes_with_sub_meshes(self):
        self.skipTest("This test is not yet implemented")
        raise NotImplementedError

    def test_2D_mesh(self):
        self.skipTest("This test is not yet implemented")
        raise NotImplementedError

    def test_2D_mesh_with_sub_meshes(self):
        self.skipTest("This test is not yet implemented")
        raise NotImplementedError

    def test_2D_mesh_in_3D_space(self):
        self.skipTest("This test is not yet implemented")
        raise NotImplementedError

    def test_2D_mesh_in_3D_space_with_sub_meshes(self):
        self.skipTest("This test is not yet implemented")
        raise NotImplementedError

    def test_3D_mesh(self):
        self.skipTest("This test is not yet implemented")
        raise NotImplementedError

    def test_3D_mesh_with_sub_meshes(self):
        self.skipTest("This test is not yet implemented")
        raise NotImplementedError

    def test_line_2N_linear_mesh(self):
        self.skipTest("This test is not yet implemented")
        raise NotImplementedError

    def test_triangle_3N_linear_mesh(self):
        self.skipTest("This test is not yet implemented")
        raise NotImplementedError

    def test_quadrilateral_4N_linear_mesh(self):
        self.skipTest("This test is not yet implemented")
        raise NotImplementedError

    def test_tetrahedra_4N_linear_mesh(self):
        def mp_check_fct(model_part):
            self.assertEqual(model_part.NumberOfNodes(), 36)
            self.assertEqual(model_part.NumberOfGeometries(), 167)

            # check how many geoms of each type
            exp_geoms = {KM.Tetrahedra3D4: 65, KM.Triangle3D3: 64, KM.Line3D2: 32, KM.Geometry: 6}
            self.assertEqual(sum(exp_geoms.values()), model_part.NumberOfGeometries())
            self.assertDictEqual(exp_geoms, get_num_geometries_by_type(model_part))

            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeLength(model_part), 3200)
            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeArea(model_part), 340000)
            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeVolume(model_part), 10000000)
            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeDomainSize(model_part), 10343200)

            for node in model_part.Nodes:
                self.assertTrue(0.0 <= node.X <= 500.0)
                self.assertTrue(0.0 <= node.X0 <= 500.0)

                self.assertTrue(0.0 <= node.Y <= 100.0)
                self.assertTrue(0.0 <= node.Y0 <= 100.0)

                self.assertTrue(0.0 <= node.Z <= 200.0)
                self.assertTrue(0.0 <= node.Z0 <= 200.0)

        self._execute_tests(Path("tetrahedral_4N"), mp_check_fct)

    def test_tetrahedra_10N_quadratic_mesh(self):
        def mp_check_fct(model_part):
            self.assertEqual(model_part.NumberOfNodes(), 168)
            self.assertEqual(model_part.NumberOfGeometries(), 176)

            # check how many geoms of each type
            exp_geoms = {KM.Tetrahedra3D10: 65, KM.Triangle3D6: 64, KM.Line3D3: 32, KM.Geometry: 15}
            self.assertEqual(sum(exp_geoms.values()), model_part.NumberOfGeometries())
            self.assertDictEqual(exp_geoms, get_num_geometries_by_type(model_part))

            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeLength(model_part), 3200)
            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeArea(model_part), 340000)
            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeVolume(model_part), 10000000)
            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeDomainSize(model_part), 10343200)

            for node in model_part.Nodes:
                self.assertTrue(0.0 <= node.X <= 500.0)
                self.assertTrue(0.0 <= node.X0 <= 500.0)

                self.assertTrue(0.0 <= node.Y <= 100.0)
                self.assertTrue(0.0 <= node.Y0 <= 100.0)

                self.assertTrue(0.0 <= node.Z <= 200.0)
                self.assertTrue(0.0 <= node.Z0 <= 200.0)

        self._execute_tests(Path("tetrahedral_10N"), mp_check_fct)

    def test_hexahedra_8N_linear_mesh(self):
        def mp_check_fct(model_part):
            self.assertEqual(model_part.NumberOfNodes(), 216)
            self.assertEqual(model_part.NumberOfGeometries(), 371)

            # check how many geoms of each type
            exp_geoms = {KM.Hexahedra3D8: 125, KM.Quadrilateral3D4: 150, KM.Line3D2: 60, KM.Geometry: 36}
            self.assertEqual(sum(exp_geoms.values()), model_part.NumberOfGeometries())
            self.assertDictEqual(exp_geoms, get_num_geometries_by_type(model_part))

            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeLength(model_part), 3200)
            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeArea(model_part), 340000)
            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeVolume(model_part), 10000000)
            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeDomainSize(model_part), 10343200)

            for node in model_part.Nodes:
                self.assertTrue(0.0 <= node.X <= 500.0)
                self.assertTrue(0.0 <= node.X0 <= 500.0)

                self.assertTrue(0.0 <= node.Y <= 100.0)
                self.assertTrue(0.0 <= node.Y0 <= 100.0)

                self.assertTrue(0.0 <= node.Z <= 200.0)
                self.assertTrue(0.0 <= node.Z0 <= 200.0)

        self._execute_tests(Path("hexahedral_8N"), mp_check_fct)

    def test_hexahedra_20N_quadratic_mesh(self):
        self.skipTest("The connectivity conversion is not yet fully implemented")

        def mp_check_fct(model_part):
            self.assertEqual(model_part.NumberOfNodes(), 36)
            self.assertEqual(model_part.NumberOfGeometries(), 36)

            # check how many geoms of each type
            exp_geoms = {KM.Hexahedra3D20: 125, KM.Quadrilateral3D8: 150, KM.Line3D3: 60, KM.Geometry: 96}
            self.assertEqual(sum(exp_geoms.values()), model_part.NumberOfGeometries())
            self.assertDictEqual(exp_geoms, get_num_geometries_by_type(model_part))

            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeLength(model_part), 3200)
            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeArea(model_part), 340000)
            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeVolume(model_part), 10000000)
            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeDomainSize(model_part), 10343200)

            for node in model_part.Nodes:
                self.assertTrue(0.0 <= node.X <= 500.0)
                self.assertTrue(0.0 <= node.X0 <= 500.0)

                self.assertTrue(0.0 <= node.Y <= 100.0)
                self.assertTrue(0.0 <= node.Y0 <= 100.0)

                self.assertTrue(0.0 <= node.Z <= 200.0)
                self.assertTrue(0.0 <= node.Z0 <= 200.0)

        self._execute_tests(Path("hexahedral_20N"), mp_check_fct)

    def test_hexahedra_27N_biquadratic_mesh(self):
        self.skipTest("The connectivity conversion is not yet fully implemented")

        def mp_check_fct(model_part):
            self.assertEqual(model_part.NumberOfNodes(), 36)
            self.assertEqual(model_part.NumberOfGeometries(), 36)

            # check how many geoms of each type
            exp_geoms = {KM.Hexahedra3D27: 125, KM.Quadrilateral3D9: 150, KM.Line3D3: 60, KM.Geometry: 121}
            self.assertEqual(sum(exp_geoms.values()), model_part.NumberOfGeometries())
            self.assertDictEqual(exp_geoms, get_num_geometries_by_type(model_part))

            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeLength(model_part), 3200)
            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeArea(model_part), 340000)
            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeVolume(model_part), 10000000)
            self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeDomainSize(model_part), 10343200)

            for node in model_part.Nodes:
                self.assertTrue(0.0 <= node.X <= 500.0)
                self.assertTrue(0.0 <= node.X0 <= 500.0)

                self.assertTrue(0.0 <= node.Y <= 100.0)
                self.assertTrue(0.0 <= node.Y0 <= 100.0)

                self.assertTrue(0.0 <= node.Z <= 200.0)
                self.assertTrue(0.0 <= node.Z0 <= 200.0)

        self._execute_tests(Path("hexahedral_27N"), mp_check_fct)


def write_vtk(model_part, name):
    # using the modeler to create elements for visualization,
    # until the vtk-output supports geometries directly
    modeler_parameters = KM.Parameters(
        """{
        "elements_list" : [{
            "model_part_name" : "read_1",
            "element_name" : "Element2D3N;Element2D4N;Element3D4N;Element3D8N;Element3D10N;Element3D20N;Element3D27N"
        }]
    }"""
    )
    modeler = KM.CreateEntitiesFromGeometriesModeler(model_part.GetModel(), modeler_parameters)
    modeler.SetupModelPart()

    vtk_parameters = KM.Parameters(
        """{
        "file_format"                  : "binary",
        "output_sub_model_parts"       : false,
        "save_output_files_in_folder"  : true
    }"""
    )

    vtk_io = KM.VtkOutput(model_part, vtk_parameters)
    vtk_io.PrintOutput(name)

    # remove elements again to avoid intereference with other checks
    KM.VariableUtils().SetFlag(KM.TO_ERASE, True, model_part.Elements)
    model_part.RemoveElementsFromAllLevels(KM.TO_ERASE)


if __name__ == "__main__":
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()
