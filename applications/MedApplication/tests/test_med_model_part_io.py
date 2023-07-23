import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.MedApplication as KratosMed

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting

from pathlib import Path



def GetMedPath(med_path, med_name="mesh.med"):
    return Path(__file__).absolute().parent / "med_files" / med_path / med_name


class TestMedModelPartIO(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()
        self.mp_read_1 = self.model.CreateModelPart("read_1")
        self.mp_read_2 = self.model.CreateModelPart("read_2")

    def _execute_tests(self, med_path, check_fct):
        med_io_read_1 = KratosMed.MedModelPartIO(GetMedPath(med_path))
        med_io_read_1.ReadModelPart(self.mp_read_1)

        with self.subTest("check_model_part"):
            check_fct(self.mp_read_1)

        with self.subTest("read_write_read"):
            med_temp_path = GetMedPath(med_path, "temp.med")
            DeleteFileIfExisting(med_temp_path) # make sure there are no leftovers from previous tests
            self.addCleanup(DeleteFileIfExisting, med_temp_path) # clean up after test

            med_io_write = KratosMed.MedModelPartIO(med_temp_path, KM.IO.WRITE)
            med_io_write.WriteModelPart(self.mp_read_1)

            med_io_read_2 = KratosMed.MedModelPartIO(med_temp_path, KM.IO.READ)
            med_io_read_2.ReadModelPart(self.mp_read_2)

            KratosMed.MedTestingUtilities.CheckModelPartsAreEqual(self.mp_read_1, self.mp_read_2)

    def test_empty_med_file(self):
        def mp_check_fct(model_part):
            self.assertEqual(model_part.NumberOfNodes(), 0)

        self._execute_tests("empty", mp_check_fct)

    def test_only_nodes(self):
        def mp_check_fct(model_part):
            self.assertEqual(model_part.NumberOfNodes(), 4)

            exp_coords = [
                (0,0,0), (0,0,1),(0,1,1),(1,1,1)
            ]

            for coords, node in zip(exp_coords, model_part.Nodes):
                self.assertAlmostEqual(node.X, coords[0])
                self.assertAlmostEqual(node.X0, coords[0])
                self.assertAlmostEqual(node.Y, coords[1])
                self.assertAlmostEqual(node.Y0, coords[1])
                self.assertAlmostEqual(node.Z, coords[2])
                self.assertAlmostEqual(node.Z0, coords[2])

        self._execute_tests("only_nodes", mp_check_fct)

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
        self.skipTest("This test is not yet implemented")
        raise NotImplementedError

    def test_hexahedra_8N_linear_mesh(self):
        self.skipTest("This test is not yet implemented")
        raise NotImplementedError


if __name__ == '__main__':
    KratosUnittest.main()
