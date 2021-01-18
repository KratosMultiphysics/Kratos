from unittest.mock import patch


try:
    import h5py
except ModuleNotFoundError:
    h5py = None


import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.HDF5Application.xdmf_utils import TryOpenH5File
from KratosMultiphysics.HDF5Application.xdmf_utils import CreateXdmfSpatialGrid
from KratosMultiphysics.HDF5Application.xdmf_utils import XdmfNodalResults
from KratosMultiphysics.HDF5Application.xdmf_utils import XdmfElementResults
from KratosMultiphysics.HDF5Application.xdmf_utils import XdmfConditionResults
from KratosMultiphysics.HDF5Application.xdmf_utils import XdmfNodalFlags
from KratosMultiphysics.HDF5Application.xdmf_utils import XdmfElementFlags
from KratosMultiphysics.HDF5Application.xdmf_utils import XdmfConditionFlags
from KratosMultiphysics.HDF5Application.xdmf_utils import XdmfResults
from KratosMultiphysics.HDF5Application.xdmf_utils import TimeLabel
from KratosMultiphysics.HDF5Application.xdmf_utils import FindMatchingFiles
from KratosMultiphysics.HDF5Application.xdmf_utils import CreateXdmfTemporalGridFromMultifile
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting


class TestTryOpenH5File(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_TryOpenH5File_NormalExecution(self):
        ok = False
        with TryOpenH5File("test.h5", "w", "core", backing_store=False) as f:
            f.create_group("grp")
            ok = True
        self.assertTrue(ok)

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_TryOpenH5File_OSError(self):
        with TryOpenH5File("test.h5", "r", "core") as f:
            self.assertTrue(f == None)

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_TryOpenH5File_KeyError(self):
        ok = False
        try:
            with TryOpenH5File("test.h5", "w", "core", backing_store=False) as f:
                f["does-not-exist"]
        except KeyError:
            ok = True
        self.assertTrue(ok)


class TestCreateXdmfSpatialGrid(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_CreateXdmfSpatialGrid(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            f.create_dataset(
                "/ModelPart/Nodes/Local/Coordinates", (15, 3), "float64")
            elem2d4n = f.create_group("/ModelPart/Xdmf/Elements/Element2D4N")
            elem2d4n.attrs["Dimension"] = 2
            elem2d4n.attrs["NumberOfNodes"] = 4
            elem2d4n.create_dataset("Connectivities", (10, 4), "int32")
            sgrid = CreateXdmfSpatialGrid(f["/ModelPart"])
        self.assertEqual(sgrid.grids[0].name, "RootModelPart.Elements.Element2D4N")
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


class TestXdmfNodalResults(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfNodalResults1(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            f.create_dataset(
                "/Results/NodalSolutionStepData/VELOCITY", (15, 3), "float64")
            results = XdmfNodalResults(f["/Results"])
        self.assertEqual(results[0].name, "VELOCITY")
        self.assertEqual(results[0].data.file_name, "kratos.h5")
        self.assertEqual(results[0].data.name, "/Results/NodalSolutionStepData/VELOCITY")
        self.assertEqual(results[0].data.dtype, "float64")
        self.assertEqual(results[0].attribute_type, "Vector")

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfNodalResults2(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            f.create_dataset(
                "/Results/NodalDataValues/DENSITY", (15,), "float64")
            results = XdmfNodalResults(f["/Results"])
        self.assertEqual(results[0].name, "DENSITY")
        self.assertEqual(results[0].data.file_name, "kratos.h5")
        self.assertEqual(results[0].data.name, "/Results/NodalDataValues/DENSITY")
        self.assertEqual(results[0].data.dtype, "float64")
        self.assertEqual(results[0].attribute_type, "Scalar")

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfNodalFlags(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            f.create_dataset(
                "/Results/NodalFlagValues/SLIP", (15,), "int32")
            results = XdmfNodalFlags(f["/Results"])
        self.assertEqual(results[0].name, "SLIP")
        self.assertEqual(results[0].data.file_name, "kratos.h5")
        self.assertEqual(results[0].data.name, "/Results/NodalFlagValues/SLIP")
        self.assertEqual(results[0].data.dtype, "int32")
        self.assertEqual(results[0].attribute_type, "Scalar")

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfNodalResults_RepeatedVariableException(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            f.create_dataset(
                "/Results/NodalSolutionStepData/DENSITY", (15,), "float64")
            f.create_dataset(
                "/Results/NodalDataValues/DENSITY", (15,), "float64")
            self.assertRaises(RuntimeError, XdmfNodalResults, f["/Results"])

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfNodalResults_NoResultsFound(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            f.create_group("/Results/NodalSolutionStepData")
            results = XdmfNodalResults(f["/Results"])
            self.assertEqual(len(results), 0)

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfNodalFlags_NoFlagsFound(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            f.create_group("/Results")
            results = XdmfNodalFlags(f["/Results"])
            self.assertEqual(len(results), 0)


class TestXdmfElementResults(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfElementResults(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            f.create_dataset(
                "/Results/ElementDataValues/DENSITY", (15,), "float64")
            results = XdmfElementResults(f["/Results"])
        self.assertEqual(results[0].name, "DENSITY")
        self.assertEqual(results[0].data.file_name, "kratos.h5")
        self.assertEqual(results[0].data.name, "/Results/ElementDataValues/DENSITY")
        self.assertEqual(results[0].data.dtype, "float64")
        self.assertEqual(results[0].attribute_type, "Scalar")

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfElementFlags(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            f.create_dataset(
                "/Results/ElementFlagValues/SLIP", (15,), "int32")
            results = XdmfElementFlags(f["/Results"])
        self.assertEqual(results[0].name, "SLIP")
        self.assertEqual(results[0].data.file_name, "kratos.h5")
        self.assertEqual(results[0].data.name, "/Results/ElementFlagValues/SLIP")
        self.assertEqual(results[0].data.dtype, "int32")
        self.assertEqual(results[0].attribute_type, "Scalar")

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfElementResults_NoResultsFound(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            f.create_group("/Results")
            results = XdmfElementResults(f["/Results"])
            self.assertEqual(len(results), 0)

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfElementFlags_NoFlagsFound(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            f.create_group("/Results")
            results = XdmfElementFlags(f["/Results"])
            self.assertEqual(len(results), 0)


class TestXdmfConditionResults(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfConditionResults(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            f.create_dataset(
                "/Results/ConditionDataValues/DENSITY", (15,), "float64")
            results = XdmfConditionResults(f["/Results"])
        self.assertEqual(results[0].name, "DENSITY")
        self.assertEqual(results[0].data.file_name, "kratos.h5")
        self.assertEqual(results[0].data.name, "/Results/ConditionDataValues/DENSITY")
        self.assertEqual(results[0].data.dtype, "float64")
        self.assertEqual(results[0].attribute_type, "Scalar")

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfConditionFlags(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            f.create_dataset(
                "/Results/ConditionFlagValues/SLIP", (15,), "int32")
            results = XdmfConditionFlags(f["/Results"])
        self.assertEqual(results[0].name, "SLIP")
        self.assertEqual(results[0].data.file_name, "kratos.h5")
        self.assertEqual(results[0].data.name, "/Results/ConditionFlagValues/SLIP")
        self.assertEqual(results[0].data.dtype, "int32")
        self.assertEqual(results[0].attribute_type, "Scalar")

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfConditionResults_NoResultsFound(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            f.create_group("/Results")
            results = XdmfConditionResults(f["/Results"])
            self.assertEqual(len(results), 0)

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfConditionFlags_NoFlagsFound(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            f.create_group("/Results")
            results = XdmfConditionFlags(f["/Results"])
            self.assertEqual(len(results), 0)


class TestXdmfResults(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_XdmfResults(self):
        with h5py.File("kratos.h5", "a", "core", backing_store=False) as f:
            f.create_group("/Results")
            results = XdmfResults(f["/Results"])
            self.assertEqual(len(results), 0)


class TestTimeLabel(KratosUnittest.TestCase):

    def test_TimeLabel(self):
        self.assertEqual(TimeLabel("kratos.h5"), "")
        self.assertEqual(TimeLabel("kratos-123.h5"), "123")
        self.assertEqual(TimeLabel("kratos-1.2.h5"), "1.2")
        self.assertEqual(TimeLabel("kratos-1.2e+00.h5"), "1.2e+00")
        self.assertEqual(TimeLabel("kratos-1.2E+00.h5"), "1.2E+00")
        self.assertEqual(TimeLabel("kratos-1.2E-01.h5"), "1.2E-01")
        self.assertEqual(TimeLabel("kratos-kratos-1.2E-01.h5"), "1.2E-01")
        self.assertEqual(TimeLabel("kratos1-kratos-1.2E-01.h5"), "1.2E-01")
        self.assertEqual(TimeLabel("kratos1-kratos-1.2E+01.h5"), "1.2E+01")


class TestFindMatchingFiles(KratosUnittest.TestCase):

    def test_FindMatchingFiles(self):
        patcher = patch("KratosMultiphysics.HDF5Application.xdmf_utils.os.listdir", autospec=True)
        listdir = patcher.start()
        listdir.return_value = ["./sim/kratos.h5",
                                "./sim/kratos-0.0000.h5", "./sim/kratos-0.2000.h5"]
        files = FindMatchingFiles("./sim/kratos")
        self.assertEqual(len(files), 3)
        self.assertTrue("./sim/kratos.h5" in files)
        self.assertTrue("./sim/kratos-0.0000.h5" in files)
        self.assertTrue("./sim/kratos-0.2000.h5" in files)
        patcher.stop()


class TestCreateXdmfTemporalGridFromMultifile(KratosUnittest.TestCase):

    def setUp(self):
        with h5py.File("kratos.h5") as f0:
            f0.create_dataset(
            "/ModelPart/Nodes/Local/Coordinates", (15, 3), "float64")
            elem2d4n = f0.create_group("/ModelPart/Xdmf/Elements/Element2D4N")
            elem2d4n.attrs["Dimension"] = 2
            elem2d4n.attrs["NumberOfNodes"] = 4
            elem2d4n.create_dataset("Connectivities", (10, 4), "int32")
        with h5py.File("kratos-1.0.h5") as f1:
            f1.create_dataset(
            "/Results/NodalSolutionStepData/VELOCITY", (15, 3), "float64")

    def tearDown(self):
        DeleteFileIfExisting("kratos.h5")
        DeleteFileIfExisting("kratos-1.0.h5")

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_CreateXdmfTemporalGridFromMultifile_OnlyOneMesh(self):
        tgrid = CreateXdmfTemporalGridFromMultifile(
            ["kratos.h5", "kratos-1.0.h5"], "/ModelPart", "/Results")
        self.assertEqual(len(tgrid.times), 2)
        self.assertEqual(len(tgrid.grids), 2)
        time0 = tgrid.times[0]
        sgrid0 = tgrid.grids[0]
        ugrid0 = sgrid0.grids[0]
        self.assertEqual(time0.time, "0.0")
        self.assertEqual(ugrid0.name, "RootModelPart.Elements.Element2D4N")
        self.assertEqual(len(ugrid0.attributes), 0)
        time1 = tgrid.times[1]
        sgrid1 = tgrid.grids[1]
        ugrid1 = sgrid1.grids[0]
        self.assertEqual(time1.time, "1.0")
        self.assertEqual(ugrid1.name, "RootModelPart.Elements.Element2D4N")
        self.assertEqual(len(ugrid1.attributes), 1)
        result1 = ugrid1.attributes[0]
        self.assertEqual(result1.name, "VELOCITY")

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_CreateXdmfTemporalGridFromMultifile_MoreThanOneMesh(self):
        with h5py.File("kratos-1.0.h5", "w") as f:
            f.create_dataset(
            "/ModelPart/Nodes/Local/Coordinates", (20, 3), "float64")
            elem2d4n = f.create_group("/ModelPart/Xdmf/Elements/Element3D3N")
            elem2d4n.attrs["Dimension"] = 3
            elem2d4n.attrs["NumberOfNodes"] = 3
            elem2d4n.create_dataset("Connectivities", (8, 3), "int32")
            f.create_dataset(
            "/Results/NodalSolutionStepData/VELOCITY", (20, 3), "float64")
        tgrid = CreateXdmfTemporalGridFromMultifile(
            ["kratos.h5", "kratos-1.0.h5"], "/ModelPart", "/Results")
        self.assertEqual(len(tgrid.times), 2)
        self.assertEqual(len(tgrid.grids), 2)
        sgrid0 = tgrid.grids[0]
        ugrid0 = sgrid0.grids[0]
        self.assertEqual(ugrid0.name, "RootModelPart.Elements.Element2D4N")
        sgrid1 = tgrid.grids[1]
        ugrid1 = sgrid1.grids[0]
        self.assertEqual(ugrid1.name, "RootModelPart.Elements.Element3D3N")

    @KratosUnittest.skipIf(h5py == None, "this test requires h5py")
    def test_CreateXdmfTemporalGridFromMultifile_XdmfNotFound(self):
        with h5py.File("kratos-1.0.h5", "w") as f:
            f.create_dataset(
            "/ModelPart/Nodes/Local/Coordinates", (20, 3), "float64")
            f.create_dataset(
            "/Results/NodalSolutionStepData/VELOCITY", (20, 3), "float64")
        tgrid = CreateXdmfTemporalGridFromMultifile(
            ["kratos.h5", "kratos-1.0.h5"], "/ModelPart", "/Results")
        self.assertEqual(len(tgrid.times), 1)
        self.assertEqual(len(tgrid.grids), 1)
        time0 = tgrid.times[0]
        self.assertEqual(time0.time, "0.0")


if __name__ == "__main__":
    KratosUnittest.main()
