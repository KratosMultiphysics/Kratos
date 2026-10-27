import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.distributed import distributed_utils
try:
    from physicsnemo.distributed.manager import DistributedManager
    have_physicsnemo = True
except ImportError:
    have_physicsnemo = False


class _Stub:
    def __init__(self, rank, size):
        self._rank, self._size = rank, size

    def Rank(self):
        return self._rank

    def Size(self):
        return self._size


class _StubManager:
    def __init__(self, rank, world_size):
        self.rank, self.world_size = rank, world_size


class TestConsistencyCheck(KratosUnittest.TestCase):
    """No process group needed: the check is pure comparison logic."""

    def test_ConsistentViewsPass(self):
        distributed_utils._CheckConsistency(_StubManager(2, 4), _Stub(2, 4))

    def test_RankMismatchRaises(self):
        with self.assertRaisesRegex(RuntimeError, "disagree"):
            distributed_utils._CheckConsistency(_StubManager(0, 4), _Stub(1, 4))

    def test_SizeMismatchRaises(self):
        with self.assertRaisesRegex(RuntimeError, "disagree"):
            distributed_utils._CheckConsistency(_StubManager(0, 2), _Stub(0, 4))


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestInitializeDistributedManager(KratosUnittest.TestCase):
    def tearDown(self):
        # Reset the singleton so repeated tests don't trip the initialized guard.
        if DistributedManager.is_initialized():
            DistributedManager.cleanup()

    def test_SerialSetup(self):
        manager = distributed_utils.InitializeDistributedManager(backend="gloo", port="29512")
        self.assertEqual(manager.rank, 0)
        self.assertEqual(manager.world_size, 1)

    def test_IdempotentWhenAlreadyInitialized(self):
        first = distributed_utils.InitializeDistributedManager(backend="gloo", port="29513")
        second = distributed_utils.InitializeDistributedManager(backend="gloo", port="29513")
        self.assertEqual(second.rank, first.rank)

    def test_AlreadyInitializedWithMismatchedCommunicatorRaises(self):
        distributed_utils.InitializeDistributedManager(backend="gloo", port="29514")
        with self.assertRaisesRegex(RuntimeError, "disagree"):
            distributed_utils.InitializeDistributedManager(data_communicator=_Stub(1, 3))


class TestGeometryNameMap(KratosUnittest.TestCase):
    """The shadow-reconstruction name maps must stay valid and complete."""

    def test_ElementMapCoversAllTessellationTypes(self):
        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import tessellation

        missing = tessellation.GetSupportedGeometryTypes() - set(
            distributed_utils._ELEMENT_NAME_BY_GEOMETRY)
        self.assertEqual(missing, set(), f"unmapped tessellation geometry types: {missing}")

    def test_EveryMappedNameConstructs(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("NameMapCheck")
        properties = model_part.CreateNewProperties(1)
        for i in range(1, 40):
            model_part.CreateNewNode(i, 0.01 * i, 0.001 * i * i, 0.1 * (i % 7))

        from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge.tessellation import _CORNER_COUNT
        entity_id = 0
        for geometry_type, name in distributed_utils._ELEMENT_NAME_BY_GEOMETRY.items():
            entity_id += 1
            node_count = int(name.rstrip("N").rsplit("D", 1)[-1])
            model_part.CreateNewElement(name, entity_id, list(range(1, node_count + 1)), properties)
        for geometry_type, name in distributed_utils._CONDITION_NAME_BY_GEOMETRY.items():
            entity_id += 1
            node_count = int(name.rstrip("N").rsplit("D", 1)[-1])
            model_part.CreateNewCondition(name, entity_id, list(range(1, node_count + 1)), properties)
        self.assertGreater(_CORNER_COUNT[next(iter(distributed_utils._ELEMENT_NAME_BY_GEOMETRY))], 0)

    def test_UnmappedConditionTypeRaises(self):
        with self.assertRaisesRegex(RuntimeError, "No generic registered condition"):
            distributed_utils._GetEntityCreationName(
                Kratos.GeometryData.KratosGeometryType.Kratos_Tetrahedra3D4, "Conditions")


class TestSerialGatherModelPartPassThrough(KratosUnittest.TestCase):
    def test_SerialPartIsReturnedUnchanged(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("Serial")
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        properties = model_part.CreateNewProperties(1)
        for i, xyz in enumerate([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]):
            model_part.CreateNewNode(i + 1, *[float(c) for c in xyz])
        model_part.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], properties)

        specs = [("PRESSURE", "node_historical"), ("PRESSURE", "element_gauss_point")]
        gathered = distributed_utils.GatherModelPartToRank0(model_part, specs)
        self.assertIs(gathered.model_part, model_part)
        self.assertIsNone(gathered.model)
        # serial pass-through leaves specs untouched (gauss translation is a
        # distributed-only concern)
        self.assertEqual(gathered.field_specs, specs)

    def test_SerialMeshGatherIsIdSorted(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("Serial")
        properties = model_part.CreateNewProperties(1)
        # deliberately created out of id order
        for node_id, xyz in [(4, (0, 0, 1)), (1, (0, 0, 0)), (3, (0, 1, 0)), (2, (1, 0, 0))]:
            model_part.CreateNewNode(node_id, *[float(c) for c in xyz])
        model_part.CreateNewElement("Element3D4N", 5, [1, 2, 3, 4], properties)

        mesh = distributed_utils.GatherMeshToRank0(model_part)
        self.assertEqual(mesh.node_ids.tolist(), [1, 2, 3, 4])
        self.assertEqual(mesh.entity_ids.tolist(), [5])
        self.assertEqual(
            Kratos.GeometryData.KratosGeometryType(int(mesh.geometry_codes[0])),
            Kratos.GeometryData.KratosGeometryType.Kratos_Tetrahedra3D4)
        self.assertEqual(mesh.connectivity[0].tolist(), [1, 2, 3, 4])


@KratosUnittest.skipUnless(have_physicsnemo, "Missing required python module: physicsnemo.")
class TestProcessGroupHelpers(KratosUnittest.TestCase):
    def tearDown(self):
        if DistributedManager.is_initialized():
            DistributedManager.cleanup()

    def test_UninitializedManagerRaises(self):
        with self.assertRaisesRegex(RuntimeError, "not initialized"):
            distributed_utils.CreateMatchedProcessGroup("g_uninitialized", 1)
        with self.assertRaisesRegex(RuntimeError, "not initialized"):
            distributed_utils.InitializeDeviceMesh((1,), ("data",))

    def test_SerialMatchedGroup(self):
        distributed_utils.InitializeDistributedManager(backend="gloo", port="29515")
        data_communicator = Kratos.ParallelEnvironment.GetDefaultDataCommunicator()
        sub = distributed_utils.CreateMatchedProcessGroup("g_serial", 1)
        self.assertIs(sub, data_communicator)  # serial: torch-side group only
        manager = DistributedManager()
        self.assertEqual(manager.group_size("g_serial"), 1)
        self.assertEqual(manager.group_rank("g_serial"), 0)

    def test_ParametersDrivenGroups(self):
        distributed_utils.InitializeDistributedManager(backend="gloo", port="29516")
        groups = distributed_utils.CreateMatchedProcessGroups(Kratos.Parameters("""{
            "process_groups" : [ { "name" : "g_parameters", "size" : 1 } ]
        }"""))
        self.assertEqual(set(groups), {"g_parameters"})
        self.assertEqual(DistributedManager().group_size("g_parameters"), 1)

    def test_DeviceMeshValidation(self):
        distributed_utils.InitializeDistributedManager(backend="gloo", port="29517")
        mesh = distributed_utils.InitializeDeviceMesh((1,), ("data",))
        self.assertIsNotNone(mesh)
        with self.assertRaisesRegex(ValueError, "implies 2 rank"):
            distributed_utils.InitializeDeviceMesh((2,), ("data",))
        with self.assertRaisesRegex(ValueError, "entries but"):
            distributed_utils.InitializeDeviceMesh((1, 1), ("data",))
        with self.assertRaisesRegex(ValueError, "may be -1"):
            distributed_utils.InitializeDeviceMesh((-1, -1), ("a", "b"))


if __name__ == '__main__':
    KratosUnittest.main()
