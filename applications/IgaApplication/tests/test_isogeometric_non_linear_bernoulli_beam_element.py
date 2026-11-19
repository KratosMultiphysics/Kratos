import numpy as np

import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import python_linear_solver_factory

# Straight cantilever along the local tangent t, with the cross section dimension WIDTH along the
# normal n and HEIGHT along v = t x n. I33 is the moment of inertia about v (deflections along n) and
# I22 the one about n (deflections along v). Torsion constant as in Bauer et al. (2016), Table 1.
E_MOD = 70000.0
NU = 0.0
G_MOD = E_MOD / (2.0 * (1.0 + NU))
WIDTH = 2.5
HEIGHT = 5.0
LENGTH = 4.0
AREA = WIDTH * HEIGHT
I33_ = HEIGHT * WIDTH ** 3 / 12.0
I22_ = WIDTH * HEIGHT ** 3 / 12.0
I_T_ = HEIGHT * WIDTH ** 3 / 3.0

DOFS = (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z, KM.ROTATION_X)
REACTIONS = (KM.REACTION_X, KM.REACTION_Y, KM.REACTION_Z, KM.REACTION_MOMENT_X)


def create_model_part(model, body_forces=False):
    model_part = model.CreateModelPart("IgaModelPart")
    for variable in (KM.DISPLACEMENT, KM.ROTATION, KM.REACTION, KM.REACTION_MOMENT, SMA.POINT_LOAD):
        model_part.AddNodalSolutionStepVariable(variable)
    if body_forces:
        model_part.AddNodalSolutionStepVariable(KM.VOLUME_ACCELERATION)
    return model_part


def create_properties(model_part, properties_id, normal, young_modulus=E_MOD, poisson_ratio=NU,
                      area=AREA, i22=I22_, i33=I33_, i_t=I_T_):
    properties = model_part.GetProperties()[properties_id]
    properties.SetValue(IGA.CROSS_AREA, area)
    properties.SetValue(IGA.I22, i22)
    properties.SetValue(IGA.I33, i33)
    properties.SetValue(IGA.I_T, i_t)
    properties.SetValue(KM.YOUNG_MODULUS, young_modulus)
    properties.SetValue(KM.POISSON_RATIO, poisson_ratio)
    properties.SetValue(KM.DENSITY, 1.0)
    properties.SetValue(IGA.LOCAL_AXIS_ORIENTATION, np.array([[0.0, *normal], [1.0, *normal]], dtype=float))
    properties.SetValue(KM.CONSTITUTIVE_LAW, IGA.BernoulliBeamElasticConstitutiveLaw())
    return properties


def create_straight_curve(model_part, start, direction, n_cp=5, degree=4, first_node_id=1):
    """Straight curve with the control points at the Greville abscissae, i.e. with a linear parametrization."""
    knots = [0.0] * degree + [i / (n_cp - degree) for i in range(1, n_cp - degree)] + [1.0] * degree
    full_knots = [0.0] + knots + [1.0]

    nodes = KM.NodesVector()
    for i in range(n_cp):
        greville_abscissa = sum(full_knots[i + 1:i + 1 + degree]) / degree
        position = np.array(start, dtype=float) + LENGTH * greville_abscissa * np.array(direction, dtype=float)
        nodes.append(model_part.CreateNewNode(first_node_id + i, *position))

    return KM.NurbsCurveGeometry3D(nodes, degree, KM.Vector(knots))


def refine_curve(model, refinements_file_name):
    modeler_parameters = KM.Parameters(r'''{
        "echo_level": 0,
        "refinements_file_name": ""
    }''')
    modeler_parameters["refinements_file_name"].SetString(refinements_file_name)
    KM.CreateModeler("RefinementModeler", model, modeler_parameters).PrepareGeometryModel()


def create_elements(model_part, curve, properties, first_element_id=1):
    quadrature_point_geometries = KM.GeometriesVector()
    curve.CreateQuadraturePointGeometries(quadrature_point_geometries, 3)
    for i in range(len(quadrature_point_geometries)):
        model_part.CreateNewElement("NonLinearBernoulliBeamElement3D", first_element_id + i, quadrature_point_geometries[i], properties)


def add_dofs(model_part):
    for dof, reaction in zip(DOFS, REACTIONS):
        KM.VariableUtils().AddDof(dof, reaction, model_part)


def clamp(nodes, axial_component):
    """Clamps an axis-aligned beam: all dofs of the first control point and the transverse displacements of the second one."""
    for dof in DOFS:
        nodes[0].Fix(dof)
    for component, dof in enumerate(DOFS[:3]):
        if component != axial_component:
            nodes[1].Fix(dof)


def add_point_load(model_part, condition_id, node, force):
    condition = model_part.CreateNewCondition("PointLoadCondition3D1N", condition_id, [node.Id], model_part.GetProperties()[100])
    condition.SetValue(SMA.POINT_LOAD, KM.Array3([float(f) for f in force]))


def add_point_moment(model_part, condition_id, node, moment):
    condition = model_part.CreateNewCondition("PointMomentCondition3D1N", condition_id, [node.Id], model_part.GetProperties()[100])
    condition.SetValue(SMA.POINT_MOMENT, KM.Array3([float(m) for m in moment]))


def solve_linear(model_part):
    """Solves with a single Newton iteration from the undeformed state, i.e. the linearized problem."""
    model_part.SetBufferSize(1)
    model_part.CloneTimeStep(1.0)

    linear_solver = python_linear_solver_factory.ConstructSolver(
        KM.Parameters(r'{"solver_type": "skyline_lu_factorization"}'))
    strategy = KM.ResidualBasedNewtonRaphsonStrategy(
        model_part,
        KM.ResidualBasedIncrementalUpdateStaticScheme(),
        linear_solver,
        KM.ResidualCriteria(1e-10, 1e-12),
        1,      # maximum iterations
        False,  # compute reactions
        False,  # reform dofs at each iteration
        False)  # move mesh
    strategy.SetEchoLevel(0)
    strategy.Check()
    strategy.Solve()


def cantilever_tip_deflection(force, moment_of_inertia):
    return force * LENGTH ** 3 / (3.0 * E_MOD * moment_of_inertia)


class BeamTestCase(KratosUnittest.TestCase):

    def assertRelativelyEqual(self, value, reference, tolerance):
        self.assertLess(abs(value - reference), tolerance * abs(reference),
                        "%.12e != %.12e (relative tolerance %.1e)" % (value, reference, tolerance))


class TestBeam4pElement(BeamTestCase):
    """Linear response of the NonLinearBernoulliBeamElement3D compared to beam theory."""

    @staticmethod
    def solve_cantilever(force, moment=(0.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0), volume_acceleration=None):
        """Cantilever along +X, clamped at x = 0 and loaded at the tip. Returns the tip node."""
        model = KM.Model()
        model_part = create_model_part(model, body_forces=volume_acceleration is not None)
        properties = create_properties(model_part, 1, normal)

        curve = create_straight_curve(model_part, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
        create_elements(model_part, curve, properties)
        IGA.ComputeBeamVectorsProcess(model_part, curve).ExecuteInitialize()
        add_dofs(model_part)

        nodes = [model_part.GetNode(i) for i in range(1, 6)]
        clamp(nodes, 0)
        if volume_acceleration is not None:
            for node in nodes:
                node.SetSolutionStepValue(KM.VOLUME_ACCELERATION, KM.Array3([float(a) for a in volume_acceleration]))
        add_point_load(model_part, 1, nodes[-1], force)
        if any(moment):
            # PointMomentCondition3D1N assembles all three moment components, but the beam only has the
            # rotation about its axis (ROTATION_X for this beam along +X). The other two are suppressed.
            for dof, reaction in ((KM.ROTATION_Y, KM.REACTION_MOMENT_Y), (KM.ROTATION_Z, KM.REACTION_MOMENT_Z)):
                nodes[-1].AddDof(dof, reaction)
                nodes[-1].Fix(dof)
            add_point_moment(model_part, 2, nodes[-1], moment)

        solve_linear(model_part)
        return nodes[-1]

    def test_constitutive_law_stress_matches_tangent(self):
        # Strains are [axial, bending n, bending v, torsion 1, torsion 2]. The stresses returned by the
        # law must be the ones of its constitutive matrix, including the two torsional components.
        model = KM.Model()
        model_part = create_model_part(model)
        properties = create_properties(model_part, 1, (0.0, 0.0, 1.0))
        geometry = KM.Line3D2(model_part.CreateNewNode(1, 0.0, 0.0, 0.0), model_part.CreateNewNode(2, LENGTH, 0.0, 0.0))

        options = KM.Flags()
        options.Set(KM.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN, True)
        options.Set(KM.ConstitutiveLaw.COMPUTE_STRESS, True)
        options.Set(KM.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)

        strain = KM.Vector([1e-3, 2e-3, -3e-3, 4e-3, -4e-3])
        stress = KM.Vector(5)
        constitutive_matrix = KM.Matrix(5, 5)
        parameters = KM.ConstitutiveLawParameters()
        parameters.SetOptions(options)
        parameters.SetStrainVector(strain)
        parameters.SetStressVector(stress)
        parameters.SetConstitutiveMatrix(constitutive_matrix)
        parameters.SetMaterialProperties(properties)
        parameters.SetElementGeometry(geometry)
        parameters.SetProcessInfo(model_part.ProcessInfo)

        constitutive_law = properties[KM.CONSTITUTIVE_LAW].Clone()
        constitutive_law.CalculateMaterialResponsePK2(parameters)

        stiffnesses = [E_MOD * AREA, E_MOD * I33_, E_MOD * I22_, 0.5 * G_MOD * I_T_, 0.5 * G_MOD * I_T_]
        expected_stress = np.array(constitutive_matrix) @ np.array(strain)
        for i, stiffness in enumerate(stiffnesses):
            self.assertRelativelyEqual(constitutive_matrix[i, i], stiffness, 1e-12)
            self.assertRelativelyEqual(stress[i], expected_stress[i], 1e-12)

    def test_axial_force(self):
        tip = self.solve_cantilever((1.0, 0.0, 0.0))
        displacement = tip.GetSolutionStepValue(KM.DISPLACEMENT)
        reference = LENGTH / (E_MOD * AREA)

        self.assertRelativelyEqual(displacement[0], reference, 1e-8)
        self.assertLess(abs(displacement[1]), 1e-8 * reference)
        self.assertLess(abs(displacement[2]), 1e-8 * reference)

    def test_bending_along_binormal(self):
        tip = self.solve_cantilever((0.0, -1.0, 0.0))
        displacement = tip.GetSolutionStepValue(KM.DISPLACEMENT)
        reference = cantilever_tip_deflection(-1.0, I22_)

        self.assertRelativelyEqual(displacement[1], reference, 1e-8)
        self.assertLess(abs(displacement[0]), 1e-8 * abs(reference))
        self.assertLess(abs(displacement[2]), 1e-8 * abs(reference))

    def test_bending_along_normal(self):
        tip = self.solve_cantilever((0.0, 0.0, -1.0))
        displacement = tip.GetSolutionStepValue(KM.DISPLACEMENT)
        reference = cantilever_tip_deflection(-1.0, I33_)

        self.assertRelativelyEqual(displacement[2], reference, 1e-8)
        self.assertLess(abs(displacement[0]), 1e-8 * abs(reference))
        self.assertLess(abs(displacement[1]), 1e-8 * abs(reference))

    def test_torsion(self):
        tip = self.solve_cantilever((0.0, 0.0, 0.0), moment=(1.0, 0.0, 0.0))
        reference = LENGTH / (G_MOD * I_T_)

        self.assertRelativelyEqual(tip.GetSolutionStepValue(KM.ROTATION_X), reference, 1e-8)
        self.assertLess(abs(tip.GetSolutionStepValue(KM.DISPLACEMENT_Y)), 1e-8 * reference)
        self.assertLess(abs(tip.GetSolutionStepValue(KM.DISPLACEMENT_Z)), 1e-8 * reference)

    def test_self_weight(self):
        # Uniform line load q = rho * A * g from the body forces: tip deflection q L^4 / (8 E I)
        gravity = 9.81
        tip = self.solve_cantilever((0.0, 0.0, 0.0), volume_acceleration=(0.0, -gravity, 0.0))
        line_load = -1.0 * AREA * gravity  # DENSITY = 1
        reference = line_load * LENGTH ** 4 / (8.0 * E_MOD * I22_)

        self.assertRelativelyEqual(tip.GetSolutionStepValue(KM.DISPLACEMENT_Y), reference, 1e-8)
        self.assertLess(abs(tip.GetSolutionStepValue(KM.DISPLACEMENT_Z)), 1e-8 * abs(reference))

    def test_rotated_cross_section(self):
        # Rotating the cross section by 90 degrees about the beam axis (normal along +Y instead of +Z)
        # swaps the bending stiffness of the two transverse directions.
        tip = self.solve_cantilever((0.0, -1.0, 0.0), normal=(0.0, 1.0, 0.0))
        self.assertRelativelyEqual(tip.GetSolutionStepValue(KM.DISPLACEMENT_Y), cantilever_tip_deflection(-1.0, I33_), 1e-8)

        tip = self.solve_cantilever((0.0, 0.0, -1.0), normal=(0.0, 1.0, 0.0))
        self.assertRelativelyEqual(tip.GetSolutionStepValue(KM.DISPLACEMENT_Z), cantilever_tip_deflection(-1.0, I22_), 1e-8)

    def test_several_beams_in_one_model_part(self):
        # Three cantilevers with different directions and cross-section orientations share one model part.
        # A single ComputeBeamVectorsProcess call must give every beam the initial frame of its own curve.
        model = KM.Model()
        model_part = create_model_part(model)

        beams = [
            # start, direction, normal, axial component, tip force, expected tip deflection
            ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0, (0.0, -1.0, 0.0), cantilever_tip_deflection(-1.0, I22_)),
            ((0.0, 10.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0, (0.0, 0.0, 1.0), cantilever_tip_deflection(1.0, I33_)),
            ((10.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), 2, (1.0, 0.0, 0.0), cantilever_tip_deflection(1.0, I33_)),
        ]

        tips = []
        for i, (start, direction, normal, axial_component, force, _) in enumerate(beams):
            properties = create_properties(model_part, i + 1, normal)
            curve = create_straight_curve(model_part, start, direction, first_node_id=10 * i + 1)
            curve.SetId(i + 1)
            model_part.AddGeometry(curve)  # keeps the parent curve of the quadrature points alive
            create_elements(model_part, curve, properties, first_element_id=100 * i + 1)
            nodes = [model_part.GetNode(10 * i + j) for j in range(1, 6)]
            tips.append((nodes, axial_component, force))

        IGA.ComputeBeamVectorsProcess(model_part).ExecuteInitialize()
        add_dofs(model_part)

        for i, (nodes, axial_component, force) in enumerate(tips):
            clamp(nodes, axial_component)
            add_point_load(model_part, i + 1, nodes[-1], force)

        solve_linear(model_part)

        for (nodes, _, force), beam in zip(tips, beams):
            loaded_component = int(np.argmax(np.abs(force)))
            deflection = nodes[-1].GetSolutionStepValue(KM.DISPLACEMENT)[loaded_component]
            self.assertRelativelyEqual(deflection, beam[5], 1e-8)

    def test_curved_quadrant_bending_torsion(self):
        # Bending-torsion interaction of a clamped quarter circle loaded out of plane at the tip,
        # Bauer et al. (2016), Section 5.2, analytical solution Eqs. (81) and (82).
        young_modulus, poisson_ratio, radius, force = 1.999e8, 0.25, 1.0, 1.0
        height, width = 0.1, 0.01
        shear_modulus = young_modulus / (2.0 * (1.0 + poisson_ratio))
        i_out_of_plane = height ** 3 * width / 12.0
        i_torsion = height * width ** 3 / 3.0

        model = KM.Model()
        model_part = create_model_part(model)
        properties = create_properties(
            model_part, 1, (0.0, 0.0, 1.0),
            young_modulus=young_modulus, poisson_ratio=poisson_ratio, area=height * width,
            i22=height * width ** 3 / 12.0, i33=i_out_of_plane, i_t=i_torsion)

        # Exact quarter circle of degree 2 from (R, 0, 0) to (0, R, 0), refined to degree 4 with 16 knot spans
        nodes = KM.NodesVector()
        for i, point in enumerate(((radius, 0.0, 0.0), (radius, radius, 0.0), (0.0, radius, 0.0)), 1):
            nodes.append(model_part.CreateNewNode(i, *point))
        curve = KM.NurbsCurveGeometry3D(nodes, 2, KM.Vector([0.0, 0.0, 1.0, 1.0]), KM.Vector([1.0, np.sqrt(0.5), 1.0]))
        curve.SetId(2)
        model_part.AddGeometry(curve)

        refine_curve(model, "modeler_tests/beam_quadrant_refinements.iga.json")

        curve = model_part.GetGeometry(2)
        create_elements(model_part, curve, properties)
        IGA.ComputeBeamVectorsProcess(model_part, curve).ExecuteInitialize()
        add_dofs(model_part)

        nodes = sorted(model_part.Nodes, key=lambda node: node.Id)
        clamp(nodes, 1)
        add_point_load(model_part, 1, nodes[-1], (0.0, 0.0, force))

        solve_linear(model_part)

        tip_displacement_reference = (force * radius ** 3 * np.pi / (4.0 * young_modulus * i_out_of_plane)
                                      + force * radius ** 3 / (shear_modulus * i_torsion) * (3.0 * np.pi / 4.0 - 2.0))
        tip_rotation_reference = (-force * radius ** 2 * np.pi / (4.0 * young_modulus * i_out_of_plane)
                                  + force * radius ** 2 / (shear_modulus * i_torsion) * (1.0 - np.pi / 4.0))

        self.assertRelativelyEqual(nodes[-1].GetSolutionStepValue(KM.DISPLACEMENT_Z), tip_displacement_reference, 1e-6)
        self.assertRelativelyEqual(nodes[-1].GetSolutionStepValue(KM.ROTATION_X), tip_rotation_reference, 1e-6)


class TestBeam4pElementNonlinear(BeamTestCase):
    """Geometrically nonlinear response of the NonLinearBernoulliBeamElement3D."""

    @staticmethod
    def build_cantilever(n_cp=5, degree=4):
        """Cantilever along +X, refined with curve_h_refinements.iga.json. Returns (model, model_part, nodes)."""
        model = KM.Model()
        model_part = create_model_part(model)
        properties = create_properties(model_part, 1, (0.0, 0.0, 1.0))

        curve = create_straight_curve(model_part, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), n_cp, degree)
        curve.SetId(2)
        model_part.AddGeometry(curve)

        refine_curve(model, "modeler_tests/curve_h_refinements.iga.json")

        curve = model_part.GetGeometry(2)
        create_elements(model_part, curve, properties)
        IGA.ComputeBeamVectorsProcess(model_part, curve).ExecuteInitialize()
        add_dofs(model_part)

        nodes = sorted(model_part.Nodes, key=lambda node: node.Id)
        clamp(nodes, 0)

        model_part.SetBufferSize(1)
        model_part.CloneTimeStep(1.0)
        for element in model_part.Elements:
            element.Initialize(model_part.ProcessInfo)
        return model, model_part, nodes

    @staticmethod
    def free_dofs(nodes):
        fixed = {4 * i + j for i, node in enumerate(nodes) for j, dof in enumerate(DOFS) if node.IsFixed(dof)}
        return [i for i in range(4 * len(nodes)) if i not in fixed]

    @staticmethod
    def element_dof_indices(element):
        # The refinement modeler numbers the control points consecutively from 1.
        return np.array([4 * (node.Id - 1) + j for node in element.GetGeometry() for j in range(4)])

    @staticmethod
    def apply_state(nodes, u):
        for i, node in enumerate(nodes):
            node.SetSolutionStepValue(KM.DISPLACEMENT, KM.Array3([float(v) for v in u[4 * i:4 * i + 3]]))
            node.SetSolutionStepValue(KM.ROTATION_X, float(u[4 * i + 3]))

    @classmethod
    def assemble(cls, model_part, nodes, u, build_level=0):
        """Global tangent and residual at the state u. BUILD_LEVEL 1/2 restricts the tangent to its material/geometric part."""
        cls.apply_state(nodes, u)
        model_part.ProcessInfo[IGA.BUILD_LEVEL] = build_level
        n_dofs = 4 * len(nodes)
        K = np.zeros((n_dofs, n_dofs))
        R = np.zeros(n_dofs)
        for element in model_part.Elements:
            lhs, rhs = KM.Matrix(), KM.Vector()
            element.CalculateLocalSystem(lhs, rhs, model_part.ProcessInfo)
            indices = cls.element_dof_indices(element)
            K[np.ix_(indices, indices)] += np.array(lhs)
            R[indices] += np.array(rhs)
        return K, R

    @classmethod
    def newton(cls, model_part, nodes, f_ext, n_steps, tolerance=1e-10, max_iterations=40):
        """Load-controlled Newton-Raphson. Returns the solution and the residual norms of the last load step."""
        free = cls.free_dofs(nodes)
        u = np.zeros(4 * len(nodes))
        for step in range(1, n_steps + 1):
            f = f_ext * step / n_steps
            scale = max(1.0, np.linalg.norm(f[free]))
            norms = []
            for _ in range(max_iterations):
                K, R = cls.assemble(model_part, nodes, u)
                residual = R + f
                norms.append(np.linalg.norm(residual[free]))
                if norms[-1] < tolerance * scale:
                    break
                u[free] += np.linalg.solve(K[np.ix_(free, free)], residual[free])
        return u, norms

    @classmethod
    def buckling_loads(cls, model_part, nodes, tip_load, reference_load):
        """Linearized buckling loads (K_mat + lambda * K_geo) after a pre-load of the given tip load."""
        f_ext = np.zeros(4 * len(nodes))
        f_ext[4 * (len(nodes) - 1):4 * (len(nodes) - 1) + 3] = tip_load
        u, _ = cls.newton(model_part, nodes, f_ext, n_steps=1)

        free = cls.free_dofs(nodes)
        K_mat, _ = cls.assemble(model_part, nodes, u, build_level=1)
        K_geo, _ = cls.assemble(model_part, nodes, u, build_level=2)

        mu = np.linalg.eigvals(-np.linalg.solve(K_mat[np.ix_(free, free)], K_geo[np.ix_(free, free)]))
        mu = mu[np.abs(mu.imag) < 1e-8 * (1.0 + np.abs(mu.real))].real
        load_factors = 1.0 / mu[np.abs(mu) > 1e-12]
        return np.sort(load_factors[load_factors > 1e-8]) * reference_load

    def test_tangent_matches_finite_differences(self):
        _, model_part, nodes = self.build_cantilever()
        n_dofs = 4 * len(nodes)

        u0 = np.zeros(n_dofs)
        for i, node in enumerate(nodes):
            x = node.X0
            u0[4 * i:4 * i + 4] = (0.004 * x, 0.010 * x ** 2, -0.006 * x ** 2, 0.003 * x)

        K, _ = self.assemble(model_part, nodes, u0)
        K_mat, _ = self.assemble(model_part, nodes, u0, build_level=1)
        K_geo, _ = self.assemble(model_part, nodes, u0, build_level=2)

        h = 1e-7
        K_fd = np.zeros_like(K)
        for j in range(n_dofs):
            u_plus = u0.copy()
            u_plus[j] += h
            u_minus = u0.copy()
            u_minus[j] -= h
            _, R_plus = self.assemble(model_part, nodes, u_plus)
            _, R_minus = self.assemble(model_part, nodes, u_minus)
            K_fd[:, j] = -(R_plus - R_minus) / (2.0 * h)

        scale = np.abs(K).max()
        self.assertLess(np.abs(K - K_fd).max() / scale, 1e-6,
                        "the tangent stiffness is not the derivative of the internal forces")
        self.assertLess(np.abs(K - (K_mat + K_geo)).max() / scale, 1e-12,
                        "the material and geometric parts do not add up to the tangent stiffness")

    def test_quadratic_newton_convergence(self):
        _, model_part, nodes = self.build_cantilever()
        f_ext = np.zeros(4 * len(nodes))
        f_ext[4 * (len(nodes) - 1) + 1] = -0.5 * E_MOD * I22_ / LENGTH ** 2

        _, norms = self.newton(model_part, nodes, f_ext, n_steps=4)
        relative_norms = np.array(norms) / norms[0]

        self.assertGreaterEqual(len(norms), 3)
        self.assertLess(relative_norms[-1], 1e-10, "Newton did not converge")
        self.assertLess(relative_norms[-2], relative_norms[-3] ** 1.8,
                        "convergence is not quadratic: %s" % " ".join("%.2e" % n for n in relative_norms))

    def test_large_deflection_elastica(self):
        # Cantilever with a tip load perpendicular to the undeformed axis, P L^2 / EI = 0.5.
        # Reference: inextensible elastica (shooting method). The linear solution would be uy = -0.6667, ux = 0.
        # The cantilever is stocky, so the axial stiffness is increased to match the inextensibility of the reference.
        ux_exact, uy_exact = -0.0636758499, -0.6485743026
        load = 0.5 * E_MOD * I22_ / LENGTH ** 2

        _, model_part, nodes = self.build_cantilever()
        model_part.GetProperties()[1].SetValue(IGA.CROSS_AREA, 1e4 * AREA)
        f_ext = np.zeros(4 * len(nodes))
        f_ext[4 * (len(nodes) - 1) + 1] = -load
        u, _ = self.newton(model_part, nodes, f_ext, n_steps=8)

        self.assertRelativelyEqual(u[4 * (len(nodes) - 1) + 1], uy_exact, 5e-4)
        self.assertRelativelyEqual(u[4 * (len(nodes) - 1)], ux_exact, 5e-4)

    def test_euler_buckling_weak_axis(self):
        _, model_part, nodes = self.build_cantilever()
        loads = self.buckling_loads(model_part, nodes, (-1000.0, 0.0, 0.0), 1000.0)
        euler = np.pi ** 2 * E_MOD * I33_ / (4.0 * LENGTH ** 2)
        self.assertRelativelyEqual(loads[0], euler, 5e-3)

    def test_euler_buckling_strong_axis(self):
        # The second mode is the first strong-axis mode, since the second weak-axis mode is at 9 times the first one.
        _, model_part, nodes = self.build_cantilever()
        loads = self.buckling_loads(model_part, nodes, (-1000.0, 0.0, 0.0), 1000.0)
        euler = np.pi ** 2 * E_MOD * I22_ / (4.0 * LENGTH ** 2)
        self.assertRelativelyEqual(loads[1], euler, 5e-3)

    def test_lateral_torsional_buckling(self):
        # Cantilever with a tip load at the centroid, Timoshenko & Gere: P_cr = 4.013 sqrt(E I_weak G I_T) / L^2
        _, model_part, nodes = self.build_cantilever()
        loads = self.buckling_loads(model_part, nodes, (0.0, -1000.0, 0.0), 1000.0)
        timoshenko = 4.013 * np.sqrt(E_MOD * I33_ * G_MOD * I_T_) / LENGTH ** 2
        self.assertRelativelyEqual(loads[0], timoshenko, 5e-2)

    def test_linear_bending_is_mesh_independent(self):
        # The exact cubic deflection is in the spline space of every mesh, since the curve is linearly parametrized.
        load = 1e-4 * E_MOD * I22_ / LENGTH ** 2
        reference = cantilever_tip_deflection(-load, I22_)
        for n_cp in (5, 7, 9, 13):
            _, model_part, nodes = self.build_cantilever(n_cp)
            f_ext = np.zeros(4 * len(nodes))
            f_ext[4 * (len(nodes) - 1) + 1] = -load
            u, _ = self.newton(model_part, nodes, f_ext, n_steps=1)
            self.assertRelativelyEqual(u[4 * (len(nodes) - 1) + 1], reference, 1e-6)


if __name__ == "__main__":
    KratosUnittest.main()
