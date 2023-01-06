import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.MedApplication as KratosMed

class TestMedModelPartIO(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def test_empty_med_file(self):
        raise NotImplementedError

    def test_only_nodes(self):
        raise NotImplementedError

    def test_nodes_with_sub_meshes(self):
        raise NotImplementedError

    def test_2D_mesh(self):
        raise NotImplementedError

    def test_2D_mesh_with_sub_meshes(self):
        raise NotImplementedError

    def test_2D_mesh_in_3D_space(self):
        raise NotImplementedError

    def test_2D_mesh_in_3D_space_with_sub_meshes(self):
        raise NotImplementedError

    def test_3D_mesh(self):
        raise NotImplementedError

    def test_3D_mesh_with_sub_meshes(self):
        raise NotImplementedError

    def test_line_2N_linear_mesh(self):
        raise NotImplementedError

    def test_triangle_3N_linear_mesh(self):
        raise NotImplementedError

    def test_quadrilateral_4N_linear_mesh(self):
        raise NotImplementedError

    def test_tetrahedra_4N_linear_mesh(self):
        raise NotImplementedError

    def test_hexahedra_8N_linear_mesh(self):
        raise NotImplementedError


if __name__ == '__main__':
    KratosUnittest.main()
