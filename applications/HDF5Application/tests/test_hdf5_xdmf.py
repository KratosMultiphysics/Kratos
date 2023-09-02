from pathlib import Path
from unittest.mock import patch


try:
    import h5py
except ModuleNotFoundError:
    h5py = None

import KratosMultiphysics as Kratos

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.HDF5Application.core.xdmf import EntityData
from KratosMultiphysics.HDF5Application.xdmf_utils import CreateXdmfSpatialGrid

class TestCreateXdmfSpatialGrid(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_CreateXdmfSpatialGrid1(self):
        self.addCleanup(DeleteFileIfExisting, str((Path(__file__).parent / f"kratos.h5").absolute()))

        with h5py.File("kratos.h5", "w") as f:
            f.create_dataset(
                "/ModelPart/Nodes/Local/Coordinates", (15, 3), "float64")
            f["/ModelPart"].attrs["__model_part_name"] = [97, 98, 99, 100]
            elem2d4n = f.create_group("/ModelPart/Xdmf/Elements/Element2D4N")
            elem2d4n.attrs["WorkingSpaceDimension"] = 2
            elem2d4n.attrs["NumberOfNodes"] = 4
            elem2d4n.create_dataset("Connectivities", (10, 4), "int32")

        sgrid = CreateXdmfSpatialGrid("kratos.h5", "/ModelPart")
        self.assertEqual(sgrid.grids[0].name, "abcd.Element2D4N")
        self.assertEqual(sgrid.grids[0].geometry.coords.file_name, "kratos.h5")
        self.assertEqual(
            sgrid.grids[0].geometry.coords.name, "/ModelPart/Nodes/Local/Coordinates")
        self.assertEqual(sgrid.grids[0].geometry.coords.dtype, "float64")
        self.assertEqual(sgrid.grids[0].geometry.coords.dimensions, (15, 3))
        self.assertEqual(sgrid.grids[0].topology.data.file_name, "kratos.h5")
        self.assertEqual(sgrid.grids[0].topology.data.name,
                         "/ModelPart/Xdmf/Elements/Element2D4N/Connectivities")
        self.assertEqual(sgrid.grids[0].topology.data.dtype, "int32")
        self.assertEqual(sgrid.grids[0].topology.data.dimensions, (10, 4))

    def test_CreateXdmfSpatialGrid2(self):
        self.addCleanup(DeleteFileIfExisting, str((Path(__file__).parent / f"kratos.h5").absolute()))

        with h5py.File("kratos.h5", "w") as f:
            f.create_dataset(
                "/ModelPart/Nodes/Local/Coordinates", (15, 3), "float64")
            f["/ModelPart"].attrs["__model_part_name"] = [97, 98, 99, 100]
            elem2d4n = f.create_group("/ModelPart/Xdmf/Elements/Element2D4N")
            elem2d4n.attrs["WorkingSpaceDimension"] = 2
            elem2d4n.attrs["NumberOfNodes"] = 4
            elem2d4n.create_dataset("Connectivities", (10, 4), "int32")

            dataset = f.create_dataset(
                "/Results/NodalSolutionStepData/VELOCITY", (15, 3), "float64")
            dataset.attrs["__mesh_location"] = [116, 101, 115, 116, 49, 46, 104, 53, 58, 47, 77, 111, 100, 101, 108, 68, 97, 116, 97]
            dataset.attrs["__data_name"] = [86, 69, 76, 79, 67, 73, 84, 89]
            dataset.attrs["__container_type"] = [78, 79, 68, 69, 83]

            dataset = f.create_dataset(
                "/Results/Nodal/VELOCITY", (15, 3), "float64")
            dataset.attrs["__mesh_location"] = [116, 101, 115, 116, 49, 46, 104, 53, 58, 47, 77, 111, 100, 101, 108, 68, 97, 116, 97]
            dataset.attrs["__data_name"] = [86, 69, 76, 79, 67, 73, 84, 89]
            dataset.attrs["__container_type"] = [78, 79, 68, 69, 83]

        sgrid = CreateXdmfSpatialGrid("kratos.h5", "/ModelPart")

        with h5py.File("kratos.h5", "r") as f:
            sgrid.AddAttribute(EntityData(f["/Results/NodalSolutionStepData/VELOCITY"]))
            with self.assertRaises(RuntimeError):
                sgrid.AddAttribute(EntityData(f["/Results/Nodal/VELOCITY"]))

    def test_CreateXdmfSpatialGrid3(self):
        self.addCleanup(DeleteFileIfExisting, str((Path(__file__).parent / f"kratos.h5").absolute()))

        with h5py.File("kratos.h5", "w") as f:
            f.create_dataset(
                "/ModelPart/Nodes/Local/Coordinates", (15, 3), "float64")
            f["/ModelPart"].attrs["__model_part_name"] = [97, 98, 99, 100]
            elem2d4n = f.create_group("/ModelPart/Xdmf/Elements/Element2D4N")
            elem2d4n.attrs["WorkingSpaceDimension"] = 2
            elem2d4n.attrs["NumberOfNodes"] = 4
            elem2d4n.create_dataset("Connectivities", (10, 4), "int32")

            dataset = f.create_dataset(
                "/Results/Test/VELOCITY", (15, 3), "float64")
            dataset.attrs["__mesh_location"] = [116, 101, 115, 116, 49, 46, 104, 53, 58, 47, 77, 111, 100, 101, 108, 68, 97, 116, 97]
            dataset.attrs["__data_name"] = [86, 69, 76, 79, 67, 73, 84, 89]
            dataset.attrs["__container_type"] = [69, 76, 69, 77, 69, 78, 84, 83]

            dataset = f.create_dataset(
                "/Results/Element/VELOCITY", (15, 3), "float64")
            dataset.attrs["__mesh_location"] = [116, 101, 115, 116, 49, 46, 104, 53, 58, 47, 77, 111, 100, 101, 108, 68, 97, 116, 97]
            dataset.attrs["__data_name"] = [86, 69, 76, 79, 67, 73, 84, 89]
            dataset.attrs["__container_type"] = [69, 76, 69, 77, 69, 78, 84, 83]

        sgrid = CreateXdmfSpatialGrid("kratos.h5", "/ModelPart")

        with h5py.File("kratos.h5", "r") as f:
            sgrid.AddAttribute(EntityData(f["/Results/Test/VELOCITY"]))
            with self.assertRaises(RuntimeError):
                sgrid.AddAttribute(EntityData(f["/Results/Element/VELOCITY"]))

class TestXdmfNodalResults(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfNodalResults1(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            dataset = f.create_dataset(
                "/Results/NodalSolutionStepData/VELOCITY", (15, 3), "float64")
            dataset.attrs["__mesh_location"] = [116, 101, 115, 116, 49, 46, 104, 53, 58, 47, 77, 111, 100, 101, 108, 68, 97, 116, 97]
            dataset.attrs["__data_name"] = [86, 69, 76, 79, 67, 73, 84, 89]
            dataset.attrs["__container_type"] = [78, 79, 68, 69, 83]
            results = EntityData(f["/Results/NodalSolutionStepData/VELOCITY"])
        self.assertEqual(results.name, "VELOCITY")
        self.assertEqual(results.container_type, Kratos.Globals.DataLocation.NodeNonHistorical)
        self.assertEqual(results.mesh_location, "test1.h5:/ModelData")
        self.assertEqual(results.data.file_name, "kratos.h5")
        self.assertEqual(results.data.name, "/Results/NodalSolutionStepData/VELOCITY")
        self.assertEqual(results.data.dtype, "float64")
        self.assertEqual(results.attribute_type, "Vector")

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfNodalResults2(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            dataset = f.create_dataset(
                "/Results/NodalDataValues/DENSITY", (15,), "float64")
            dataset.attrs["__mesh_location"] = [116, 101, 115, 116, 49, 46, 104, 53, 58, 47, 77, 111, 100, 101, 108, 68, 97, 116, 97]
            dataset.attrs["__data_name"] = [68, 69, 78, 83, 73, 84, 89]
            dataset.attrs["__container_type"] = [78, 79, 68, 69, 83]
            results = EntityData(f["/Results/NodalDataValues/DENSITY"])
        self.assertEqual(results.name, "DENSITY")
        self.assertEqual(results.container_type, Kratos.Globals.DataLocation.NodeNonHistorical)
        self.assertEqual(results.mesh_location, "test1.h5:/ModelData")
        self.assertEqual(results.data.file_name, "kratos.h5")
        self.assertEqual(results.data.name, "/Results/NodalDataValues/DENSITY")
        self.assertEqual(results.data.dtype, "float64")
        self.assertEqual(results.attribute_type, "Scalar")

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfNodalFlags(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            dataset = f.create_dataset(
                "/Results/NodalFlagValues/SLIP", (15,), "int32")
            dataset.attrs["__mesh_location"] = [116, 101, 115, 116, 49, 46, 104, 53, 58, 47, 77, 111, 100, 101, 108, 68, 97, 116, 97]
            dataset.attrs["__data_name"] = [83, 76, 73, 80]
            dataset.attrs["__container_type"] = [78, 79, 68, 69, 83]
            results = EntityData(f["/Results/NodalFlagValues/SLIP"])
        self.assertEqual(results.name, "SLIP")
        self.assertEqual(results.container_type, Kratos.Globals.DataLocation.NodeNonHistorical)
        self.assertEqual(results.mesh_location, "test1.h5:/ModelData")
        self.assertEqual(results.data.file_name, "kratos.h5")
        self.assertEqual(results.data.name, "/Results/NodalFlagValues/SLIP")
        self.assertEqual(results.data.dtype, "int32")
        self.assertEqual(results.attribute_type, "Scalar")


class TestXdmfElementResults(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfElementResults(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            dataset = f.create_dataset(
                "/Results/ElementDataValues/DENSITY", (15,), "float64")
            dataset.attrs["__mesh_location"] = [116, 101, 115, 116, 49, 46, 104, 53, 58, 47, 77, 111, 100, 101, 108, 68, 97, 116, 97]
            dataset.attrs["__data_name"] = [68, 69, 78, 83, 73, 84, 89]
            dataset.attrs["__container_type"] = [69, 76, 69, 77, 69, 78, 84, 83]
            results = EntityData(f["/Results/ElementDataValues/DENSITY"])
        self.assertEqual(results.name, "DENSITY")
        self.assertEqual(results.container_type, Kratos.Globals.DataLocation.Element)
        self.assertEqual(results.mesh_location, "test1.h5:/ModelData")
        self.assertEqual(results.data.file_name, "kratos.h5")
        self.assertEqual(results.data.name, "/Results/ElementDataValues/DENSITY")
        self.assertEqual(results.data.dtype, "float64")
        self.assertEqual(results.attribute_type, "Scalar")

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfElementFlags(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            dataset = f.create_dataset(
                "/Results/ElementFlagValues/SLIP", (15,), "int32")
            dataset.attrs["__mesh_location"] = [116, 101, 115, 116, 49, 46, 104, 53, 58, 47, 77, 111, 100, 101, 108, 68, 97, 116, 97]
            dataset.attrs["__data_name"] = [83, 76, 73, 80]
            dataset.attrs["__container_type"] = [69, 76, 69, 77, 69, 78, 84, 83]
            results = EntityData(f["/Results/ElementFlagValues/SLIP"])
        self.assertEqual(results.name, "SLIP")
        self.assertEqual(results.container_type, Kratos.Globals.DataLocation.Element)
        self.assertEqual(results.mesh_location, "test1.h5:/ModelData")
        self.assertEqual(results.data.file_name, "kratos.h5")
        self.assertEqual(results.data.name, "/Results/ElementFlagValues/SLIP")
        self.assertEqual(results.data.dtype, "int32")
        self.assertEqual(results.attribute_type, "Scalar")

class TestXdmfConditionResults(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfConditionResults(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            dataset = f.create_dataset(
                "/Results/ConditionDataValues/DENSITY", (15,), "float64")
            dataset.attrs["__mesh_location"] = [116, 101, 115, 116, 49, 46, 104, 53, 58, 47, 77, 111, 100, 101, 108, 68, 97, 116, 97]
            dataset.attrs["__data_name"] = [68, 69, 78, 83, 73, 84, 89]
            dataset.attrs["__container_type"] = [67, 79, 78, 68, 73, 84, 73, 79, 78, 83]
            results = EntityData(f["/Results/ConditionDataValues/DENSITY"])
        self.assertEqual(results.name, "DENSITY")
        self.assertEqual(results.container_type, Kratos.Globals.DataLocation.Condition)
        self.assertEqual(results.mesh_location, "test1.h5:/ModelData")
        self.assertEqual(results.data.file_name, "kratos.h5")
        self.assertEqual(results.data.name, "/Results/ConditionDataValues/DENSITY")
        self.assertEqual(results.data.dtype, "float64")
        self.assertEqual(results.attribute_type, "Scalar")

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfConditionFlags(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            dataset = f.create_dataset(
                "/Results/ConditionFlagValues/SLIP", (15,), "int32")
            dataset.attrs["__mesh_location"] = [116, 101, 115, 116, 49, 46, 104, 53, 58, 47, 77, 111, 100, 101, 108, 68, 97, 116, 97]
            dataset.attrs["__data_name"] = [83, 76, 73, 80]
            dataset.attrs["__container_type"] = [67, 79, 78, 68, 73, 84, 73, 79, 78, 83]
            results = EntityData(f["/Results/ConditionFlagValues/SLIP"])
        self.assertEqual(results.name, "SLIP")
        self.assertEqual(results.container_type, Kratos.Globals.DataLocation.Condition)
        self.assertEqual(results.mesh_location, "test1.h5:/ModelData")
        self.assertEqual(results.data.file_name, "kratos.h5")
        self.assertEqual(results.data.name, "/Results/ConditionFlagValues/SLIP")
        self.assertEqual(results.data.dtype, "int32")
        self.assertEqual(results.attribute_type, "Scalar")

if __name__ == "__main__":
    KratosUnittest.main()
