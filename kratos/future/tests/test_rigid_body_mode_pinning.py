import numpy as np

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

# This file lives in the "future" test set, which is also run by core-only builds,
# but the system is assembled with structural elements, so the application is optional
try:
    import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
    use_structural_mechanics = True
except ImportError:
    use_structural_mechanics = False


@KratosUnittest.skipUnless(use_structural_mechanics, "StructuralMechanicsApplication is required to assemble the test system")
class TestRigidBodyModePinning(KratosUnittest.TestCase):

    def _create_model_part(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2

        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        props = mp.GetProperties()[1]
        props.SetValue(KratosMultiphysics.YOUNG_MODULUS, 210e9)
        props.SetValue(KratosMultiphysics.POISSON_RATIO, 0.3)
        props.SetValue(KratosMultiphysics.THICKNESS, 1.0)
        props.SetValue(KratosMultiphysics.DENSITY, 1.0)
        props.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw())

        # unit square split in two triangles sharing the face 1-3
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        mp.CreateNewNode(3, 1.0, 1.0, 0.0)
        mp.CreateNewNode(4, 0.0, 1.0, 0.0)
        # separate unit square discretised by one quadrilateral
        mp.CreateNewNode(5, 2.0, 0.0, 0.0)
        mp.CreateNewNode(6, 3.0, 0.0, 0.0)
        mp.CreateNewNode(7, 3.0, 1.0, 0.0)
        mp.CreateNewNode(8, 2.0, 1.0, 0.0)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)

        mp.CreateNewElement("SmallDisplacementElement2D3N", 1, [1, 2, 3], props)
        mp.CreateNewElement("SmallDisplacementElement2D3N", 2, [1, 3, 4], props)
        mp.CreateNewElement("SmallDisplacementElement2D4N", 3, [5, 6, 7, 8], props)

        return mp

    def _prescribe_state(self, mp):
        # arbitrary displacement field, including a rigid translation of each component
        for node in mp.Nodes:
            x = node.X0
            y = node.Y0
            shift = 1.0e-3 if node.Id <= 4 else -2.0e-3
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, shift + 1.0e-4 * x + 3.0e-5 * y)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, -5.0e-5 * x + 2.0e-4 * y)

    def _rigid_body_vectors(self, mp, node_ids, n_dofs):
        tx = np.zeros(n_dofs)
        ty = np.zeros(n_dofs)
        rotation = np.zeros(n_dofs)
        for node in mp.Nodes:
            if node.Id not in node_ids:
                continue
            x0, y0 = node.X0, node.Y0
            dx_eqid = node.GetDof(KratosMultiphysics.DISPLACEMENT_X).EquationId
            dy_eqid = node.GetDof(KratosMultiphysics.DISPLACEMENT_Y).EquationId
            tx[dx_eqid] = 1.0
            ty[dy_eqid] = 1.0
            rotation[dx_eqid] = -y0
            rotation[dy_eqid] = x0
        return [tx, ty, rotation]

    def _expected_solution(self, mp, n_dofs):
        """The exact solution the pinning solver has to return, in closed form.

        The assembled rhs is b = -K u (the element routine provides the internal forces),
        with u the prescribed state, so the solutions of K x = b are exactly -u + z with z
        a rigid body mode of the free components, and the pinning solver must return the
        minimal-norm element of that set: -u with the rigid content of each connected
        component projected out (no reference solver needed).
        """
        u = np.zeros(n_dofs)
        for node in mp.Nodes:
            u[node.GetDof(KratosMultiphysics.DISPLACEMENT_X).EquationId] = -node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0)
            u[node.GetDof(KratosMultiphysics.DISPLACEMENT_Y).EquationId] = -node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)

        for node_ids in ((1, 2, 3, 4), (5, 6, 7, 8)):
            modes = np.array(self._rigid_body_vectors(mp, node_ids, n_dofs))  # (3 x n)
            coeffs = np.linalg.solve(modes @ modes.T, modes @ u)
            u -= coeffs @ modes
        return u

    def _build_system(self, mp, process_info, scheme, dof_set):
        """Assembles the global system directly into a CsrMatrix."""
        scheme.Initialize(mp)
        scheme.InitializeElements(mp)

        self._prescribe_state(mp)

        # the local stiffness of a structural element is fully populated,
        # so every (i, j) pair of the element dofs goes into the pattern
        equation_ids = [element.EquationIdVector(process_info) for element in mp.Elements]
        n_dofs = len(dof_set)

        graph = KratosMultiphysics.SparseContiguousRowGraph(n_dofs)
        for ids in equation_ids:
            graph.AddEntries(ids)
        graph.Finalize()

        # the assembly is additive, so A must be zero before assembling
        A = KratosMultiphysics.CsrMatrix(graph)
        A.SetValue(0.0)
        b = KratosMultiphysics.SystemVector(n_dofs)
        b.SetValue(0.0)

        K = KratosMultiphysics.Matrix()
        f = KratosMultiphysics.Vector()
        for element, ids in zip(mp.Elements, equation_ids):
            element.CalculateLocalSystem(K, f, process_info)
            A.Assemble(K, ids)
            b.Assemble(f, ids)

        return A, b, n_dofs

    def test_rigid_body_mode_pinning_free_floating_bodies(self):
        mp = self._create_model_part()
        process_info = mp.ProcessInfo

        # dof set and equation ids (the dof set is needed as additional data of the linear system)
        dof_set = KratosMultiphysics.DofsArrayType()
        KratosMultiphysics.BlockBuildDofArrayUtility.SetUpDofArray(mp, dof_set)
        KratosMultiphysics.DofArrayUtilities.SetDofEquationIds(dof_set)

        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

        A, b, n_dofs = self._build_system(mp, process_info, scheme, dof_set)
        x_pinned = KratosMultiphysics.SystemVector(n_dofs)

        # keep a copy of the original system: the solver modifies A and b in place
        A_ref = A.copy()
        b_ref = b.copy()

        # right hand side b = K u_ref belongs to the range of K, so the system is consistent
        self.assertGreater(max(abs(b_ref[i]) for i in range(n_dofs)), 0.0)

        linear_solver = KratosMultiphysics.Future.RigidBodyModePinningSolver(
            KratosMultiphysics.Parameters("""{
                "solver_type": "rigid_body_mode_pinning",
                "inner_solver_settings": {
                    "solver_type": "skyline_lu_factorization"
                },
                "tolerance": 1.0e-6,
                "echo_level": 1
            }"""))
        linear_system = KratosMultiphysics.Future.LinearSystem(A, b, x_pinned, "rigid_body_pinning")
        linear_system.SetAdditionalData(mp, dof_set)

        linear_solver.Initialize(linear_system)
        linear_solver.PrepareAdditionalData(linear_system)
        linear_solver.InitializeSolutionStep(linear_system)
        self.assertTrue(linear_solver.PerformSolutionStep(linear_system))

        # the pinning must have actually modified the matrix (zeroed pinned off-diagonal entries)
        self.assertNotEqual(list(A.value_data()), list(A_ref.value_data()))

        # residual of the original (unpinned) system must vanish: the pinned solve must
        # still satisfy the original equilibrium equations
        residual = np.array(A_ref @ x_pinned - b_ref)
        rhs_norm = np.linalg.norm(np.array(b_ref))
        self.assertLess(np.linalg.norm(residual) / rhs_norm, 1.0e-8)

        # the solution must be the minimal-norm one, which here is known in closed form:
        # with b = -K u the solution set is -u + null(K), so the correct answer is -u with
        # the rigid content of each free component projected out. Checking against it
        # verifies at once that the equilibrium is satisfied and that no arbitrary rigid
        # body content polluted the result (a singular system has no unique solution, so
        # comparing against another solver's output would have no well-defined meaning).
        x_expected = self._expected_solution(mp, n_dofs)
        expected_norm = np.linalg.norm(x_expected)

        # the test problem must have a non-trivial solution
        self.assertGreater(expected_norm, 1.0e-6)

        # the solver must reproduce it (which, together with the non-triviality bound
        # above, also guarantees the computed solution is not the trivial zero one)
        self.assertLess(np.linalg.norm(np.array(x_pinned) - x_expected), 1.0e-8 * expected_norm)

        # the result must not depend on the starting point: repeat the solve from a
        # polluted initial guess (the converged solution plus rigid content added to
        # both free components, which also contaminates the pinned dofs)
        polluted = np.array(x_pinned)
        for node_ids, scale in (((1, 2, 3, 4), 5.0e-3), ((5, 6, 7, 8), -3.0e-3)):
            tx, ty, rz = self._rigid_body_vectors(mp, node_ids, n_dofs)
            polluted += scale * (tx + ty + rz)

        x_pinned_restart = KratosMultiphysics.SystemVector(n_dofs)
        for i in range(n_dofs):
            x_pinned_restart[i] = polluted[i]

        A_restart = A_ref.copy()
        b_restart = b_ref.copy()
        system_restart = KratosMultiphysics.Future.LinearSystem(A_restart, b_restart, x_pinned_restart, "rigid_body_pinning_restart")
        system_restart.SetAdditionalData(mp, dof_set)
        linear_solver.InitializeSolutionStep(system_restart)
        self.assertTrue(linear_solver.PerformSolutionStep(system_restart))
        self.assertLess(np.linalg.norm(np.array(x_pinned_restart) - x_expected), 1.0e-8 * expected_norm)

    def _convert_to_current_interface(self, A_csr, b_system, n_dofs):
        """Builds the system arrays of the current interface from the CSR system.

        The CompressedMatrix is constructed from the CSR arrays (row pointers,
        column indices and values) of the future-built CsrMatrix.
        """
        A = KratosMultiphysics.CompressedMatrix(
            A_csr.index1_data(), A_csr.index2_data(), A_csr.value_data())
        b = KratosMultiphysics.Vector(b_system)
        return A, b

    def test_rigid_body_mode_pinning_current_interface(self):
        """Same system as test_rigid_body_mode_pinning_free_floating_bodies, but solved through the current (old) linear solver interface."""
        mp = self._create_model_part()
        process_info = mp.ProcessInfo

        # dof set and equation ids (the dof set is needed as additional data of the linear system)
        dof_set = KratosMultiphysics.DofsArrayType()
        KratosMultiphysics.BlockBuildDofArrayUtility.SetUpDofArray(mp, dof_set)
        KratosMultiphysics.DofArrayUtilities.SetDofEquationIds(dof_set)

        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

        # the system is built with the future mechanism and then converted to
        # the current interface (CompressedMatrix + Vector) via its CSR arrays
        A_csr, b_csr, n_dofs = self._build_system(mp, process_info, scheme, dof_set)
        A, b = self._convert_to_current_interface(A_csr, b_csr, n_dofs)
        x_pinned = KratosMultiphysics.Vector(n_dofs)
        KratosMultiphysics.SparseSpace().SetToZeroVector(x_pinned)

        # keep a copy of the original system: the solver modifies A and b in place
        A_ref = KratosMultiphysics.CompressedMatrix(A)
        b_ref = KratosMultiphysics.Vector(b)

        # right hand side b = K u_ref belongs to the range of K, so the system is consistent
        self.assertGreater(max(abs(b_ref[i]) for i in range(n_dofs)), 0.0)

        linear_solver = KratosMultiphysics.LinearSolverFactory().Create(
            KratosMultiphysics.Parameters("""{
                "solver_type": "rigid_body_mode_pinning",
                "inner_solver_settings" : {
                    "solver_type": "skyline_lu_factorization"
                },
                "tolerance": 1.0e-6,
                "echo_level": 1
            }"""))
        linear_solver.Initialize(A, x_pinned, b)
        linear_solver.ProvideAdditionalData(A, x_pinned, b, dof_set, mp)
        linear_solver.InitializeSolutionStep(A, x_pinned, b)
        self.assertTrue(linear_solver.PerformSolutionStep(A, x_pinned, b))

        # the pinning must have actually modified the matrix (zeroed pinned off-diagonal entries)
        self.assertNotEqual(list(A.value_data()), list(A_ref.value_data()))

        # residual of the original (unpinned) system must vanish
        residual = KratosMultiphysics.Vector(n_dofs)
        KratosMultiphysics.SparseSpace().SetToZeroVector(residual)
        KratosMultiphysics.SparseSpace().Mult(A_ref, x_pinned, residual)
        KratosMultiphysics.SparseSpace().UnaliasedAdd(residual, -1.0, b_ref) # residual -= b
        residual_norm = np.linalg.norm(np.array(residual))

        rhs_norm = np.linalg.norm(np.array(b_ref))
        self.assertLess(residual_norm / rhs_norm, 1.0e-8)

        # the solution must be the minimal-norm one, known in closed form (see
        # _expected_solution): this checks correctness and the absence of arbitrary rigid
        # body content in one shot, without comparing against another solver
        x_expected = self._expected_solution(mp, n_dofs)
        expected_norm = np.linalg.norm(x_expected)

        # the test problem must have a non-trivial solution
        self.assertGreater(expected_norm, 1.0e-6)

        # the solver must reproduce it (which, together with the non-triviality bound
        # above, also guarantees the computed solution is not the trivial zero one)
        self.assertLess(np.linalg.norm(np.array(x_pinned) - x_expected), 1.0e-8 * expected_norm)

        # the result must not depend on the starting point: repeat the solve from a
        # polluted initial guess (the converged solution plus rigid content added to
        # both free components, which also contaminates the pinned dofs)
        polluted = np.array(x_pinned)
        for node_ids, scale in (((1, 2, 3, 4), 5.0e-3), ((5, 6, 7, 8), -3.0e-3)):
            tx, ty, rz = self._rigid_body_vectors(mp, node_ids, n_dofs)
            polluted += scale * (tx + ty + rz)

        x_pinned_restart = KratosMultiphysics.Vector(n_dofs)
        for i in range(n_dofs):
            x_pinned_restart[i] = polluted[i]

        A_restart = KratosMultiphysics.CompressedMatrix(A_ref)
        b_restart = KratosMultiphysics.Vector(b_ref)
        linear_solver.InitializeSolutionStep(A_restart, x_pinned_restart, b_restart)
        self.assertTrue(linear_solver.PerformSolutionStep(A_restart, x_pinned_restart, b_restart))
        self.assertLess(np.linalg.norm(np.array(x_pinned_restart) - x_expected), 1.0e-8 * expected_norm)

    def test_parameters_constructor(self):
        # the default constructor builds the default inner solver (skyline LU)
        self.assertIsNotNone(KratosMultiphysics.Future.RigidBodyModePinningSolver())
        self.assertIsNotNone(KratosMultiphysics.Future.RigidBodyModePinningSolver(
            KratosMultiphysics.Parameters('{"inner_solver_settings": {"solver_type": "eigen_sparse_lu"}, "tolerance": 1.0e-4}')))
        # an unknown inner solver type must raise an error
        with self.assertRaisesRegex(RuntimeError, "unknown"):
            KratosMultiphysics.Future.RigidBodyModePinningSolver(
                KratosMultiphysics.Parameters('{"inner_solver_settings": {"solver_type": "not_a_solver"}}'))


if __name__ == '__main__':
    KratosUnittest.main()
