import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.MedApplication as KratosMed

from pathlib import Path
from typing import Any, Dict, List, Set


def GetMedPath(med_file_name: Path) -> Path:
    return Path(__file__).absolute().parent / "med_files" / med_file_name.with_suffix(".med")


class MedModelPartIOTestCase(KratosUnittest.TestCase):
    def _basic_checks(self, model_part):
        # check no elements or conditions are created
        self.assertEqual(model_part.NumberOfElements(), 0)
        self.assertEqual(model_part.NumberOfConditions(), 0)

        self.assertGreaterEqual(KratosMed.MedTestingUtilities.ComputeLength(model_part), 0.0)
        self.assertGreaterEqual(KratosMed.MedTestingUtilities.ComputeArea(model_part), 0.0)
        self.assertGreaterEqual(KratosMed.MedTestingUtilities.ComputeVolume(model_part), 0.0)
        self.assertGreaterEqual(KratosMed.MedTestingUtilities.ComputeDomainSize(model_part), 0.0)

        # check increasing node Ids
        for i, node in enumerate(model_part.Nodes):
            self.assertEqual(node.Id, i + 1)

        # check geometries have correct Ids
        # Note: Geometries are not ordered, thus cannot check like nodes
        self.assertEqual(set(get_geometry_ids(model_part)), set(range(1, model_part.NumberOfGeometries() + 1)))

        # check that the entities are unique in the ModelParts
        self._check_unique_nodes(model_part)
        self._check_unique_geometries(model_part)

        # check each ModelPart has (at least) the nodes of its geometries
        self._check_nodes_geometries(model_part)

    def _check_unique_nodes(self, model_part):
        node_ids: List[int] = get_node_ids(model_part)
        node_ids_unique: Set[int] = set(node_ids)
        self.assertEqual(len(node_ids), len(node_ids_unique), msg=f"Name of ModelPart: {model_part.FullName()}")

        for smp in model_part.SubModelParts:
            self._check_unique_nodes(smp)

    def _check_unique_geometries(self, model_part):
        geom_ids: List[int] = get_geometry_ids(model_part)
        geom_ids_unique: Set[int] = set(geom_ids)
        self.assertEqual(len(geom_ids), len(geom_ids_unique), msg=f"Name of ModelPart: {model_part.FullName()}")

        for smp in model_part.SubModelParts:
            self._check_unique_geometries(smp)

    def _check_nodes_geometries(self, model_part):
        geom_node_ids: List[int] = get_geom_node_ids(model_part)
        geom_node_ids_unique: Set[int] = set(geom_node_ids)

        node_ids: List[int] = get_node_ids(model_part)

        # there can be nodes that do not belong to geometries, thus geom_node_ids is subset
        self.assertTrue(geom_node_ids_unique.issubset(node_ids), msg=f"Name of ModelPart: {model_part.FullName()}")

        for smp in model_part.SubModelParts:
            self._check_nodes_geometries(smp)


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


def get_geometry_ids(model_part: KM.ModelPart) -> List[int]:
    return [geom.Id for geom in model_part.Geometries]


def get_geom_node_ids(model_part: KM.ModelPart) -> List[int]:
    node_ids: List[int] = []
    for geom in model_part.Geometries:
        for node in geom:
            node_ids.append(node.Id)
    return node_ids
