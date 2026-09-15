"""Nonlinear displacement mapping from straight and curved IGA beams to FEM boxes."""

import math
from pathlib import Path

import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication  # Register mappers.
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable


class _IgaBeamMapperBoxBase(KratosUnittest.TestCase):
    curved = False

    @classmethod
    def setUpClass(cls):
        if not CheckIfApplicationsAvailable("IgaApplication"):
            raise KratosUnittest.SkipTest("IgaApplication is required")
        import KratosMultiphysics.IgaApplication as IGA
        cls.iga = IGA

    def setUp(self):
        self.model = KM.Model()
        self.beam = self.model.CreateModelPart("beam")
        self.box = self.model.CreateModelPart("box")
        self.beam.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.beam.AddNodalSolutionStepVariable(KM.ROTATION)
        self.box.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.box.AddNodalSolutionStepVariable(KM.MESH_DISPLACEMENT)
        if "force" in self._testMethodName:
            self.skipTestIfApplicationsNotAvailable("StructuralMechanicsApplication")
            import KratosMultiphysics.StructuralMechanicsApplication as SMA
            self.sma = SMA
            self.beam.AddNodalSolutionStepVariable(SMA.POINT_LOAD)
            self.beam.AddNodalSolutionStepVariable(SMA.POINT_MOMENT)
            self.box.AddNodalSolutionStepVariable(KM.FORCE)
        self.length = 10.0
        points = KM.NodesVector()
        if self.curved:
            self.radius = 10.0
            coordinates = [(self.radius, 0.0, 0.0),
                           (self.radius, self.radius, 0.0), (0.0, self.radius, 0.0)]
        else:
            coordinates = [(i * self.length / 2, 0.0, 0.0) for i in range(3)]
        for i, point in enumerate(coordinates):
            points.append(self.beam.CreateNewNode(i + 1, *point))
        knots = KM.Vector([0.0, 0.0, 1.0, 1.0])
        if self.curved:
            self.curve = KM.NurbsCurveGeometry3D(
                points, 2, knots, KM.Vector([1.0, math.sqrt(0.5), 1.0]))
        else:
            self.curve = KM.NurbsCurveGeometry3D(points, 2, knots)
        properties = self.beam.CreateNewProperties(1)
        initial_tangent = [0.0, 1.0, 0.0] if self.curved else [1.0, 0.0, 0.0]
        properties.SetValue(self.iga.T_0, KM.Vector(initial_tangent))
        properties.SetValue(self.iga.N_0, KM.Vector([0.0, 0.0, 1.0]))
        orientation = KM.Matrix(2, 4)
        for i in range(2):
            for j in range(4):
                orientation[i, j] = 0.0
        orientation[1, 0] = 1.0
        orientation[0, 3] = orientation[1, 3] = 1.0
        properties.SetValue(self.iga.LOCAL_AXIS_ORIENTATION, orientation)
        geometries = KM.GeometriesVector()
        self.curve.CreateQuadraturePointGeometries(geometries, 3)
        for i in range(len(geometries)):
            self.beam.CreateNewElement("IsogeometricBeamElement", i + 1, geometries[i], properties)

        # Refine all six faces, sharing nodes at edges and corners.
        nx, ny, nz = 100, 4, 4
        node_ids = {}
        for i in range(nx + 1):
            for j in range(ny + 1):
                for k in range(nz + 1):
                    if i in (0, nx) or j in (0, ny) or k in (0, nz):
                        node_id = len(node_ids) + 1
                        node_ids[i, j, k] = node_id
                        x, y, z = self.length * i / nx, -0.5 + j / ny, -0.3 + 0.6 * k / nz
                        if self.curved:
                            angle = math.pi * i / (2 * nx)
                            x, y = (self.radius - y) * math.cos(angle), (self.radius - y) * math.sin(angle)
                        self.box.CreateNewNode(node_id, x, y, z)
        surface_properties = self.box.CreateNewProperties(1)
        faces = []
        for i in (0, nx):
            for j in range(ny):
                for k in range(nz):
                    face = [(i, j, k), (i, j + 1, k), (i, j + 1, k + 1), (i, j, k + 1)]
                    faces.append(face[::-1] if i == 0 else face)
        for j in (0, ny):
            for i in range(nx):
                for k in range(nz):
                    face = [(i, j, k), (i + 1, j, k), (i + 1, j, k + 1), (i, j, k + 1)]
                    faces.append(face if j == 0 else face[::-1])
        for k in (0, nz):
            for i in range(nx):
                for j in range(ny):
                    face = [(i, j, k), (i + 1, j, k), (i + 1, j + 1, k), (i, j + 1, k)]
                    faces.append(face[::-1] if k == 0 else face)
        for i, face in enumerate(faces):
            self.box.CreateNewCondition("SurfaceCondition3D4N", i + 1,
                                        [node_ids[index] for index in face], surface_properties)
        self.mapper = KM.MapperFactory.CreateMapper(
            self.beam, self.box, KM.Parameters('{"mapper_type":"iga_beam_mapper"}'))
        self.initialize_vtk_output()
        self.write_vtk_output()

    def initialize_vtk_output(self):
        from KratosMultiphysics.vtk_output_process import VtkOutputProcess
        from KratosMultiphysics.IgaApplication.iga_vtk_output_process import IgaVTKOutputProcess

        output_directory = Path(__file__).resolve().parent / "iga_beam_mapper_vtk" / self._testMethodName
        output_directory.mkdir(parents=True, exist_ok=True)
        self.vtk_curve = KM.BrepCurve(self.curve)
        self.vtk_curve.SetId(1)
        self.beam.AddGeometry(self.vtk_curve)

        surface_settings = KM.Parameters('''{
            "model_part_name": "box",
            "file_format": "ascii",
            "output_precision": 12,
            "output_path": "",
            "entity_type": "condition",
            "write_deformed_configuration": false,
            "nodal_solution_step_data_variables": ["DISPLACEMENT"],
            "output_control_type": "step",
            "output_interval": 1
        }''')
        surface_settings["output_path"].SetString(str(output_directory / "surface"))
        beam_settings = KM.Parameters('''{
            "model_part_name": "beam",
            "output_file_name": "",
            "brep_curve_ids": [1],
            "output_refinement_curve": [2],
            "nodal_solution_step_data_variables": ["DISPLACEMENT", "ROTATION"],
            "output_control_type": "step",
            "output_interval": 1
        }''')
        beam_settings["output_file_name"].SetString(str(output_directory / "beam"))
        if self.curved:
            beam_settings["output_refinement_curve"][0].SetInt(99)
        self.output_processes = [VtkOutputProcess(self.model, surface_settings),
                                 IgaVTKOutputProcess(self.model, beam_settings)]
        self.output_step = 0
        for process in self.output_processes:
            process.ExecuteInitialize()
            process.ExecuteBeforeSolutionLoop()

    def write_vtk_output(self):
        # Both writers store reference geometry plus displacement. The first
        # output caches the undeformed IGA geometry before any mesh-motion tests.
        for model_part in (self.beam, self.box):
            model_part.ProcessInfo[KM.STEP] = self.output_step
            model_part.ProcessInfo[KM.TIME] = float(self.output_step)
        for process in self.output_processes:
            process.ExecuteInitializeSolutionStep()
            process.ExecuteBeforeOutputStep()
            process.PrintOutput()
            process.ExecuteAfterOutputStep()
            process.ExecuteFinalizeSolutionStep()
        self.output_step += 1

    def tearDown(self):
        for process in self.output_processes:
            process.ExecuteFinalize()

    def check_displacements(self, expected_position):
        self.mapper.Map(KM.DISPLACEMENT, KM.DISPLACEMENT)
        for node in self.box.Nodes:
            initial = (node.X0, node.Y0, node.Z0)
            expected = expected_position(*initial)
            displacement = node.GetSolutionStepValue(KM.DISPLACEMENT)
            for i in range(3):
                self.assertAlmostEqual(displacement[i], expected[i] - initial[i], delta=1e-8)
        self.write_vtk_output()


    def initialize_load_conditions(self):
        for node in self.beam.Nodes:
            for variable in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z,
                             KM.ROTATION_X, KM.ROTATION_Y, KM.ROTATION_Z):
                node.AddDof(variable)
            node.Fix(KM.ROTATION_Y)
            node.Fix(KM.ROTATION_Z)
            self.beam.CreateNewCondition("PointLoadCondition3D1N", 2 * node.Id - 1,
                                         [node.Id], self.beam.GetProperties()[1])
            self.beam.CreateNewCondition("PointMomentCondition3D1N", 2 * node.Id,
                                         [node.Id], self.beam.GetProperties()[1])

    def check_beam_transfer_operator(self):
        self.skipTestIfApplicationsNotAvailable("CoSimulationApplication")
        from KratosMultiphysics.CoSimulationApplication.factories import data_transfer_operator_factory
        from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData

        self.initialize_load_conditions()
        settings = KM.Parameters('''{
            "type": "kratos_beam_mapping",
            "solver_name_beam": "structure",
            "solver_name_surface": "fluid",
            "model_part_name_beam": "beam",
            "model_part_name_surface": "box",
            "mapper_settings": {"mapper_type": "iga_beam_mapper"}
        }''')

        def interface(model_part, variable, solver):
            model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3
            parameters = KM.Parameters('{"model_part_name":"", "variable_name":"", "dimension":3}')
            parameters["model_part_name"].SetString(model_part.Name)
            parameters["variable_name"].SetString(variable.Name())
            return CouplingInterfaceData(parameters, self.model, solver_name=solver)

        beam_displacement = interface(self.beam, KM.DISPLACEMENT, "structure")
        surface_displacement = interface(self.box, KM.MESH_DISPLACEMENT, "fluid")
        beam_load = interface(self.beam, self.sma.POINT_LOAD, "structure")
        surface_force = interface(self.box, KM.FORCE, "fluid")
        for node in self.beam.Nodes:
            node.SetSolutionStepValue(KM.DISPLACEMENT, KM.Vector([0.1 * node.Id, 0.2, 0.05 * node.Id**2]))
            node.SetSolutionStepValue(KM.ROTATION_X, 0.4 * node.Id)
        for node in self.box.Nodes:
            node.SetSolutionStepValue(KM.FORCE, KM.Vector([0.1, math.sin(node.Id), 0.3]))
        self.mapper.Map(KM.DISPLACEMENT, KM.DISPLACEMENT)
        expected_displacement = [list(n.GetSolutionStepValue(KM.DISPLACEMENT)) for n in self.box.Nodes]
        self.mapper.InverseMap(self.sma.POINT_LOAD, KM.FORCE)
        expected_load = [list(n.GetSolutionStepValue(self.sma.POINT_LOAD)) for n in self.beam.Nodes]
        expected_moment = [list(n.GetSolutionStepValue(self.sma.POINT_MOMENT)) for n in self.beam.Nodes]
        for inverse_first in (False, True):
            operator = data_transfer_operator_factory.CreateDataTransferOperator(
                settings.Clone(), KM.Testing.GetDefaultDataCommunicator())
            transfers = [(beam_displacement, surface_displacement), (surface_force, beam_load)]
            if inverse_first:
                transfers.reverse()
            for node in self.box.Nodes:
                node.SetSolutionStepValue(KM.DISPLACEMENT, KM.Vector([0.0, 0.0, 0.0]))
            for node in self.beam.Nodes:
                node.SetSolutionStepValue(self.sma.POINT_LOAD, KM.Vector([0.0, 0.0, 0.0]))
                node.SetSolutionStepValue(self.sma.POINT_MOMENT, KM.Vector([0.0, 0.0, 0.0]))
            for source, target in transfers:
                operator.TransferData(source, target, KM.Parameters('[]'))
            for node, value in zip(self.box.Nodes, expected_displacement):
                self.assertVectorAlmostEqual(node.GetSolutionStepValue(KM.MESH_DISPLACEMENT), KM.Vector(value))
            for node, force, moment in zip(self.beam.Nodes, expected_load, expected_moment):
                self.assertVectorAlmostEqual(node.GetSolutionStepValue(self.sma.POINT_LOAD), KM.Vector(force))
                self.assertVectorAlmostEqual(node.GetSolutionStepValue(self.sma.POINT_MOMENT), KM.Vector(moment))
            operator.TransferData(surface_force, beam_load, KM.Parameters('["swap_sign"]'))
            for node, force, moment in zip(self.beam.Nodes, expected_load, expected_moment):
                self.assertVectorAlmostEqual(node.GetSolutionStepValue(self.sma.POINT_LOAD), KM.Vector([-v for v in force]))
                self.assertVectorAlmostEqual(node.GetSolutionStepValue(self.sma.POINT_MOMENT), KM.Vector([-v for v in moment]))
        with self.assertRaisesRegex(RuntimeError, "invalid secondary rotation"):
            operator.beam_mapper.Map(KM.DISPLACEMENT, KM.DISPLACEMENT, KM.DISPLACEMENT)
        with self.assertRaisesRegex(RuntimeError, "invalid secondary moment"):
            operator.beam_mapper.InverseMap(self.sma.POINT_LOAD, KM.FORCE, KM.FORCE)

    def check_force_tangent(self):
        self.initialize_load_conditions()
        for node in self.beam.Nodes:
            node.SetSolutionStepValue(KM.DISPLACEMENT, KM.Vector(
                [0.1 * node.Id, -0.07 * node.Id**2, 0.13 * node.Id**2]))
            node.SetSolutionStepValue(KM.ROTATION_X, 0.7 * node.Id)
        for node in self.box.Nodes:
            node.SetSolutionStepValue(KM.FORCE, KM.Vector(
                [math.sin(node.Id), math.cos(0.3 * node.Id), 0.2]))
        self.mapper.InverseMap(self.sma.POINT_LOAD, KM.FORCE)
        expected = []
        for node in self.beam.Nodes:
            expected.append(list(node.GetSolutionStepValue(self.sma.POINT_LOAD))
                            + [node.GetSolutionStepValue(self.sma.POINT_MOMENT)[0]])
        # Each DOF direction checks H^T f against the derivative of f dot u.
        epsilon = 1e-6
        for a, node in enumerate(self.beam.Nodes):
            for j, variable in enumerate((KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y,
                                           KM.DISPLACEMENT_Z, KM.ROTATION_X)):
                value = node.GetSolutionStepValue(variable)
                work = []
                for sign in (1, -1):
                    node.SetSolutionStepValue(variable, value + sign * epsilon)
                    self.mapper.Map(KM.DISPLACEMENT, KM.DISPLACEMENT)
                    work.append(sum(sum(n.GetSolutionStepValue(KM.FORCE)[k]
                                        * n.GetSolutionStepValue(KM.DISPLACEMENT)[k]
                                        for k in range(3)) for n in self.box.Nodes))
                node.SetSolutionStepValue(variable, value)
                self.assertAlmostEqual((work[0] - work[1]) / (2 * epsilon),
                                       expected[a][j], delta=2e-5)
        # Repeated mapping overwrites; native conditions assemble exactly these loads.
        self.mapper.InverseMap(self.sma.POINT_LOAD, KM.FORCE)
        for a, node in enumerate(self.beam.Nodes):
            for offset, variable in ((-1, self.sma.POINT_LOAD), (0, self.sma.POINT_MOMENT)):
                condition = self.beam.GetCondition(2 * node.Id + offset)
                rhs = KM.Vector()
                condition.CalculateRightHandSide(rhs, self.beam.ProcessInfo)
                self.assertVectorAlmostEqual(rhs, node.GetSolutionStepValue(variable))
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(self.sma.POINT_LOAD),
                                         KM.Vector(expected[a][:3]))
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(self.sma.POINT_MOMENT),
                                         KM.Vector([expected[a][3], 0.0, 0.0]))
        for k in range(3):
            self.assertAlmostEqual(sum(n.GetSolutionStepValue(self.sma.POINT_LOAD)[k] for n in self.beam.Nodes),
                                   sum(n.GetSolutionStepValue(KM.FORCE)[k] for n in self.box.Nodes), delta=1e-9)
        with self.assertRaisesRegex(RuntimeError, "only SWAP_SIGN"):
            self.mapper.InverseMap(self.sma.POINT_LOAD, KM.FORCE, KM.Mapper.ADD_VALUES)
        last = list(self.box.Nodes)[-1]
        last.SetSolutionStepValue(KM.FORCE, KM.Vector([float("nan"), 0.0, 0.0]))
        with self.assertRaisesRegex(RuntimeError, "non-finite surface force"):
            self.mapper.InverseMap(self.sma.POINT_LOAD, KM.FORCE)
        for a, node in enumerate(self.beam.Nodes):
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(self.sma.POINT_LOAD),
                                         KM.Vector(expected[a][:3]))
        for node in self.box.Nodes:
            node.SetSolutionStepValue(KM.FORCE, KM.Vector([0.0, 0.0, 0.0]))
        self.mapper.InverseMap(self.sma.POINT_LOAD, KM.FORCE)
        for node in self.beam.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(self.sma.POINT_LOAD), KM.Vector([0.0, 0.0, 0.0]))
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(self.sma.POINT_MOMENT), KM.Vector([0.0, 0.0, 0.0]))


class TestIgaBeamMapperBox(_IgaBeamMapperBoxBase):
    def test_endpoint_force_mapping(self):
        # Exercise every generalized-load derivative with offsets beyond both ends.
        next_id = max(node.Id for node in self.box.Nodes) + 1
        self.box.CreateNewNode(next_id, -1.0, 0.4, 0.2)
        self.box.CreateNewNode(next_id + 1, 12.0, -0.3, 0.1)
        self.mapper = KM.MapperFactory.CreateMapper(
            self.beam, self.box, KM.Parameters('{"mapper_type":"iga_beam_mapper"}'))
        self.check_force_tangent()

    def test_force_transfer_operator(self):
        self.check_beam_transfer_operator()

    def test_force_mapping(self):
        self.check_force_tangent()

    def test_pure_torque_force(self):
        self.initialize_load_conditions()
        section = [node for node in self.box.Nodes if abs(node.X0 - 5.0) < 1e-12]
        first = next(node for node in section if abs(node.Y0 - 0.5) < 1e-12 and abs(node.Z0) < 1e-12)
        second = next(node for node in section if abs(node.Y0 + 0.5) < 1e-12 and abs(node.Z0) < 1e-12)
        first.SetSolutionStepValue(KM.FORCE, KM.Vector([0.0, 0.0, 2.0]))
        second.SetSolutionStepValue(KM.FORCE, KM.Vector([0.0, 0.0, -2.0]))
        self.mapper.InverseMap(self.sma.POINT_LOAD, KM.FORCE)
        for node, weight in zip(self.beam.Nodes, (0.25, 0.5, 0.25)):
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(self.sma.POINT_LOAD), KM.Vector([0.0, 0.0, 0.0]))
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(self.sma.POINT_MOMENT),
                                         KM.Vector([2.0 * weight, 0.0, 0.0]))

    def test_eccentric_point_force(self):
        self.initialize_load_conditions()
        node = self.box.GetNode(1)
        node.SetSolutionStepValue(KM.FORCE, KM.Vector([0.0, 0.0, 2.0]))
        self.mapper.InverseMap(self.sma.POINT_LOAD, KM.FORCE)
        self.assertAlmostEqual(sum(n.GetSolutionStepValue(self.sma.POINT_MOMENT)[0]
                                   for n in self.beam.Nodes), 2 * node.Y0)
        self.assertAlmostEqual(sum(n.GetSolutionStepValue(self.sma.POINT_LOAD)[2]
                                   for n in self.beam.Nodes), 2.0)

    def test_zero_and_translation(self):
        self.assertEqual(self.box.NumberOfNodes(), 1634)
        self.assertEqual(self.box.NumberOfConditions(), 1632)
        self.check_displacements(lambda x, y, z: (x, y, z))
        for node in self.beam.Nodes:
            node.SetSolutionStepValue(KM.DISPLACEMENT, KM.Vector([0.7, -0.2, 0.4]))
        expected = lambda x, y, z: (x + 0.7, y - 0.2, z + 0.4)
        self.check_displacements(expected)
        # Ignore mesh motion and overwrite old output, rather than accumulating it.
        for node in self.box.Nodes:
            node.X += 100.0
            node.SetSolutionStepValue(KM.DISPLACEMENT, KM.Vector([9.0, 9.0, 9.0]))
        for node in self.beam.Nodes:
            node.Y += 100.0
        self.check_displacements(expected)
        self.check_displacements(expected)

    def test_twist(self):
        for node in self.beam.Nodes:
            node.SetSolutionStepValue(KM.ROTATION_X, 2 * math.pi * node.X0 / self.length)

        def expected(x, y, z):
            angle = 2 * math.pi * x / self.length
            return x, y * math.cos(angle) - z * math.sin(angle), y * math.sin(angle) + z * math.cos(angle)

        self.check_displacements(expected)

    def test_bending(self):
        tip = 0.5
        self.beam.GetNode(3).SetSolutionStepValue(KM.DISPLACEMENT_Y, tip)

        def expected(x, y, z):
            angle = math.atan(2 * tip * x / self.length**2)
            return x - y * math.sin(angle), tip * (x / self.length)**2 + y * math.cos(angle), z

        self.check_displacements(expected)

    def test_finite_rigid_rotation(self):
        bend = 2.0  # More than 90 degrees.
        twist = 4 * math.pi
        c, s = math.cos(bend), math.sin(bend)
        for node in self.beam.Nodes:
            node.SetSolutionStepValue(KM.DISPLACEMENT,
                                      KM.Vector([(c - 1) * node.X0 + 0.4, s * node.X0 - 0.3, 0.2]))
            node.SetSolutionStepValue(KM.ROTATION_X, twist)

        def expected(x, y, z):
            rotated_y = y * math.cos(twist) - z * math.sin(twist)
            rotated_z = y * math.sin(twist) + z * math.cos(twist)
            return c * x - s * rotated_y + 0.4, s * x + c * rotated_y - 0.3, rotated_z + 0.2

        self.check_displacements(expected)

    def test_invalid_states_and_flags(self):
        with self.assertRaisesRegex(RuntimeError, "without flags"):
            self.mapper.Map(KM.DISPLACEMENT, KM.DISPLACEMENT, KM.Mapper.ADD_VALUES)
        for node in self.beam.Nodes:
            node.SetSolutionStepValue(KM.DISPLACEMENT_X, -2 * node.X0)
        with self.assertRaisesRegex(RuntimeError, "opposite reference and current tangents"):
            self.mapper.Map(KM.DISPLACEMENT, KM.DISPLACEMENT)
        for node in self.beam.Nodes:
            node.SetSolutionStepValue(KM.DISPLACEMENT_X, -node.X0)
        with self.assertRaisesRegex(RuntimeError, "invalid current beam state"):
            self.mapper.Map(KM.DISPLACEMENT, KM.DISPLACEMENT)
        for node in self.box.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(KM.DISPLACEMENT), KM.Vector([0.0, 0.0, 0.0]))


class TestIgaBeamMapperCurvedBox(_IgaBeamMapperBoxBase):
    """Exact rational quarter-circle and an annular box of fixed section size."""

    curved = True

    def test_curved_force_transfer_operator(self):
        self.check_beam_transfer_operator()

    def test_curved_force_mapping(self):
        self.check_force_tangent()

    def test_curved_zero_and_translation(self):
        self.assertEqual(self.box.NumberOfNodes(), 1634)
        self.assertEqual(self.box.NumberOfConditions(), 1632)
        self.check_displacements(lambda x, y, z: (x, y, z))
        for node in self.beam.Nodes:
            node.SetSolutionStepValue(KM.DISPLACEMENT, KM.Vector([0.7, -0.2, 0.4]))
        expected = lambda x, y, z: (x + 0.7, y - 0.2, z + 0.4)
        self.check_displacements(expected)
        # Reference attachments must remain valid after mesh coordinates change.
        for node in self.beam.Nodes:
            node.X += 100.0
        for node in self.box.Nodes:
            node.Z -= 100.0
        self.check_displacements(expected)

    def test_curved_in_plane_rigid_rotation(self):
        # A rotation about Z is the bending transport for this planar arc.
        # It therefore requires no additional scalar twist.
        angle = 2.0
        c, s = math.cos(angle), math.sin(angle)

        def expected(x, y, z):
            return c * x - s * y + 0.4, s * x + c * y - 0.3, z + 0.2

        for node in self.beam.Nodes:
            initial = (node.X0, node.Y0, node.Z0)
            current = expected(*initial)
            node.SetSolutionStepValue(KM.DISPLACEMENT,
                                      KM.Vector([current[i] - initial[i] for i in range(3)]))
        self.check_displacements(expected)

    def test_curved_twist(self):
        # Rationally interpolate root/middle/tip twists 0, pi, 2*pi.
        for i, node in enumerate(self.beam.Nodes):
            node.SetSolutionStepValue(KM.ROTATION_X, i * math.pi)

        def expected(x, y, z):
            radial_distance = math.hypot(x, y)
            c, s = x / radial_distance, y / radial_distance
            offset = radial_distance - self.radius
            # For this rational circle: N0+N1=cos(theta), N1+N2=sin(theta).
            # Thus pi*N1+2*pi*N2 = pi*(1+sin(theta)-cos(theta)).
            twist = math.pi * (1 + s - c)
            current_offset = offset * math.cos(twist) + z * math.sin(twist)
            current_z = z * math.cos(twist) - offset * math.sin(twist)
            return (self.radius + current_offset) * c, (self.radius + current_offset) * s, current_z

        self.check_displacements(expected)

    def test_curved_radius_expansion(self):
        # Increase the centerline radius by 20%, keeping the section dimensions.
        factor = 1.2
        for node in self.beam.Nodes:
            node.SetSolutionStepValue(KM.DISPLACEMENT,
                                      KM.Vector([(factor - 1) * node.X0, (factor - 1) * node.Y0, 0.0]))

        def expected(x, y, z):
            radial_distance = math.hypot(x, y)
            increment = (factor - 1) * self.radius
            return x + increment * x / radial_distance, y + increment * y / radial_distance, z

        self.check_displacements(expected)


if __name__ == "__main__":
    KratosUnittest.main()
