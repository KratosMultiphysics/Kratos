import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.MedApplication as KratosMed

from pathlib import Path


def GetMedPath(med_path, med_name="mesh.med"):
    return Path(__file__).absolute().parent / "med_files" / med_path / med_name


class TestMedModelPartIOReadSubModelPart(KratosUnittest.TestCase):

    def test_cube_with_groups(self):
        model = KM.Model()
        model_part = model.CreateModelPart("test")
        KratosMed.MedModelPartIO(GetMedPath("cube_with_groups", "cube.med")).ReadModelPart(model_part)

        self.assertEqual(model_part.NumberOfNodes(), 130)
        self.assertEqual(model_part.NumberOfGeometries(), 644)
        self.assertEqual(model_part.NumberOfElements(), 0)
        self.assertEqual(model_part.NumberOfConditions(), 0)

        self.assertEqual(model_part.NumberOfSubModelParts(), 2)
        self.assertTrue(model_part.HasSubModelPart("interface"))
        self.assertTrue(model_part.HasSubModelPart("interface_nodes"))

        smp_interface = model_part.GetSubModelPart("interface")
        smp_interface_nodes = model_part.GetSubModelPart("interface_nodes")

        self.assertEqual(smp_interface.NumberOfNodes(), smp_interface_nodes.NumberOfNodes())
        self.assertEqual(smp_interface.NumberOfNodes(), 27)

        self.assertEqual(smp_interface.NumberOfGeometries(), 36)
        self.assertEqual(smp_interface_nodes.NumberOfGeometries(), 0)

        # check if smp "interface" contains exactly the nodes of its geometries
        # check coords of the nodes in the smps (x=200, 0<=y<=200, 0<=z<=200)
        # check that the correct geoms (3D triangles) are in the smp
        # check how many geoms of each type:
        #   48 geometries of type Line3D2
        #   86 geometries of type Tetrahedra3D4
        #   210 geometries of type Triangle3D3


if __name__ == '__main__':
    KratosUnittest.main()
