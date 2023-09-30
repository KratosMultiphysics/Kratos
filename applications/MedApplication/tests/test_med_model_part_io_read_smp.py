import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.MedApplication as KratosMed

from itertools import chain
from pathlib import Path
from typing import Any, Dict, List, Set


def GetMedPath(med_path, med_name="mesh.med"):
    return Path(__file__).absolute().parent / "med_files" / med_path / med_name


class TestMedModelPartIOReadSubModelPart(KratosUnittest.TestCase):

    def test_cube_with_groups(self):
        model: KM.Model = KM.Model()
        model_part: KM.ModelPart = model.CreateModelPart("test")
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

        # make sure the nodes are unique
        smp_interface_node_ids: List[int] = get_node_ids(smp_interface)
        smp_interface_node_ids_unique: Set[int] = set(smp_interface_node_ids)
        self.assertEqual(len(smp_interface_node_ids), len(smp_interface_node_ids_unique))

        smp_interface_nodes_node_ids: List[int] = get_node_ids(smp_interface_nodes)
        smp_interface_nodes_node_ids_unique: Set[int] = set(smp_interface_nodes_node_ids)
        self.assertEqual(len(smp_interface_nodes_node_ids), len(smp_interface_nodes_node_ids_unique))

        # make sure that the smp has the nodes of its geometries
        self.assertEqual(smp_interface_node_ids_unique, set(get_geom_node_ids(smp_interface)))

        # check if smps have same nodes
        self.assertEqual(smp_interface_node_ids_unique, smp_interface_nodes_node_ids_unique)

        # make sure the geometries are unique
        smp_interface_geom_ids: List[int] = [geom.Id for geom in model_part.Geometries]
        smp_interface_geom_ids_unique: Set[int] = set(smp_interface_geom_ids)
        self.assertEqual(len(smp_interface_geom_ids), len(smp_interface_geom_ids_unique))

        # check increasing node Ids
        for i, node in enumerate(model_part.Nodes):
            self.assertEqual(node.Id, i+1)

        # check geometries have correct Ids
        # Note: Geometries are not ordered, thus cannot check like nodes
        self.assertEqual(smp_interface_geom_ids_unique, set(range(1, model_part.NumberOfGeometries()+1)))

        # check node coordinates
        for node in chain(smp_interface.Nodes, smp_interface_nodes.Nodes):
            self.assertAlmostEqual(node.X, 200.0)
            self.assertAlmostEqual(node.X0, 200.0)

            self.assertTrue(0.0 <= node.Y <= 200.0)
            self.assertTrue(0.0 <= node.Y0 <= 200.0)

            self.assertTrue(0.0 <= node.Z <= 200.0)
            self.assertTrue(0.0 <= node.Z0 <= 200.0)

        # check that the correct geoms (3D triangles) are in the smp
        for geom in smp_interface.Geometries:
            self.assertIsInstance(geom, KM.Triangle3D3)

        # check how many geoms of each type
        exp_geoms: Dict[Any, int] = {
            KM.Tetrahedra3D4: 386,
            KM.Triangle3D3: 210,
            KM.Line3D2: 48
        }
        self.assertEqual(sum(exp_geoms.values()), 644)
        self.assertDictEqual(exp_geoms, get_num_geometries_by_type(model_part))


def get_num_geometries_by_type(model_part: KM.ModelPart) -> Dict[Any, int]:
    geoms_by_type: Dict[Any, int] = {}
    for geom in model_part.Geometries:
        type_geom = type(geom)
        if type_geom not in geoms_by_type:
            geoms_by_type[type_geom] = 0
        geoms_by_type[type_geom] += 1
    return geoms_by_type


def get_node_ids(model_part: KM.ModelPart) -> List[int]:
    return [node.Id for node in model_part.Nodes]


def get_geom_node_ids(model_part: KM.ModelPart) -> List[int]:
    node_ids: List[int] = []
    for geom in model_part.Geometries:
        for node in geom:
            node_ids.append(node.Id)
    return node_ids

if __name__ == '__main__':
    KratosUnittest.main()
