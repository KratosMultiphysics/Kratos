import math
import numpy as np

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest


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
        dx = KratosMultiphysics.SystemVector(n_dofs)

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
        linear_system = KratosMultiphysics.Future.LinearSystem(A, b, dx, "rigid_body_pinning")
        linear_system.SetAdditionalData(mp, dof_set)

        linear_solver.Initialize(linear_system)
        linear_solver.PrepareAdditionalData(linear_system)
        linear_solver.InitializeSolutionStep(linear_system)
        self.assertTrue(linear_solver.PerformSolutionStep(linear_system))

        # the pinning must have actually modified the matrix (zeroed pinned off-diagonal entries)
        self.assertNotEqual(list(A.value_data()), list(A_ref.value_data()))

        for i in range(n_dofs):
            self.assertTrue(abs(dx[i]) < 1.0e10)

        # residual of the original (unpinned) system must vanish
        residual = np.array(A_ref @ dx - b_ref)
        residual_norm = np.linalg.norm(residual) #math.sqrt(sum(residual[i] * residual[i] for i in range(n_dofs)))

        rhs_norm = np.linalg.norm(np.array(b_ref)) 
        self.assertLess(residual_norm / rhs_norm, 1.0e-8)

        # A is singular (each free component contributes 3 rigid-body modes), so a plain
        # LU solve does not detect the singularity and returns a solution polluted by the
        # arbitrary rigid-body content, whose norm is much larger than the one of the
        # pinned system. The pinning solver removes the rigid-body modes, so the solution
        # it produces must be significantly smaller.
        pinned_norm = np.linalg.norm(dx) 

        # the solution must be independent of the dofs that were pinned: its rigid body
        # content (the 3 rigid body vectors of each free component) must be filtered out
        dx_array = np.array(dx)
        for node_ids in ((1, 2, 3, 4), (5, 6, 7, 8)):
            for mode in self._rigid_body_vectors(mp, node_ids, n_dofs):
                self.assertLess(abs(float(dx_array @ mode)), 1.0e-9 * pinned_norm)

        dx_unpinned = KratosMultiphysics.SystemVector(n_dofs)
        unpinned_solver = KratosMultiphysics.Future.SkylineLUFactorizationSolver()
        unpinned_system = KratosMultiphysics.Future.LinearSystem(
            A_ref, b_ref, dx_unpinned, "unpinned_reference")
        unpinned_solver.Initialize(unpinned_system)
        self.assertTrue(unpinned_solver.PerformSolutionStep(unpinned_system))
        
        unpinned_norm = np.linalg.norm(np.array(dx_unpinned))
        self.assertGreater(unpinned_norm, 2.0 * pinned_norm)

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
        x = KratosMultiphysics.Vector(n_dofs)
        KratosMultiphysics.SparseSpace().SetToZeroVector(x)

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
        linear_solver.Initialize(A, x, b)
        linear_solver.ProvideAdditionalData(A, x, b, dof_set, mp)
        linear_solver.InitializeSolutionStep(A, x, b)
        self.assertTrue(linear_solver.PerformSolutionStep(A, x, b))

        # the pinning must have actually modified the matrix (zeroed pinned off-diagonal entries)
        self.assertNotEqual(list(A.value_data()), list(A_ref.value_data()))

        for i in range(n_dofs):
            self.assertTrue(abs(x[i]) < 1.0e10)

        # residual of the original (unpinned) system must vanish
        residual = KratosMultiphysics.Vector(n_dofs)
        KratosMultiphysics.SparseSpace().SetToZeroVector(residual)
        KratosMultiphysics.SparseSpace().Mult(A_ref, x, residual)
        KratosMultiphysics.SparseSpace().UnaliasedAdd(residual, -1.0, b_ref) # residual -= b
        residual_norm = np.linalg.norm(np.array(residual))

        rhs_norm = np.linalg.norm(np.array(b_ref))
        self.assertLess(residual_norm / rhs_norm, 1.0e-8)

        # A is singular (each free component contributes 3 rigid-body modes), so a plain
        # LU solve does not detect the singularity and returns a solution polluted by the
        # arbitrary rigid-body content, whose norm is much larger than the one of the
        # pinned system. The pinning solver removes the rigid-body modes, so the solution
        # it produces must be significantly smaller.
        pinned_norm = np.linalg.norm(np.array(x))

        # the solution must be independent of the dofs that were pinned: its rigid body
        # content (the 3 rigid body vectors of each free component) must be filtered out
        dx_array = np.array(x)
        for node_ids in ((1, 2, 3, 4), (5, 6, 7, 8)):
            for mode in self._rigid_body_vectors(mp, node_ids, n_dofs):
                self.assertLess(abs(float(dx_array @ mode)), 1.0e-9 * pinned_norm)

        # reference: the plain LU solve of the unpinned system produces a solution of a
        # much larger norm, since it contains the arbitrary rigid body content
        dx_unpinned = KratosMultiphysics.Vector(n_dofs)
        KratosMultiphysics.SparseSpace().SetToZeroVector(dx_unpinned)
        unpinned_solver = KratosMultiphysics.LinearSolverFactory().Create(
            KratosMultiphysics.Parameters("""{"solver_type": "skyline_lu_factorization"}"""))
        A_unpinned = KratosMultiphysics.CompressedMatrix(A_ref)
        b_unpinned = KratosMultiphysics.Vector(b_ref)
        unpinned_solver.Solve(A_unpinned, dx_unpinned, b_unpinned)
        unpinned_norm = np.linalg.norm(np.array(dx_unpinned))
        self.assertGreater(unpinned_norm, 2.0 * pinned_norm)

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
