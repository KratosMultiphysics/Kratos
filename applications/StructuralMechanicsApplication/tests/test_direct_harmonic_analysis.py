import numpy as np

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication


class DirectHarmonicAnalysisTests(KratosUnittest.TestCase):

    def _add_variables(self, mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD_IMAGINARY)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.DISPLACEMENT_IMAGINARY)

    def _add_dofs(self, mp):
        vu = KratosMultiphysics.VariableUtils()
        vu.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
        vu.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)
        vu.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, mp)

    def _create_3dof_geometry(self, current_model, stiffness, mass):
        mp = current_model.CreateModelPart("mdof_direct")
        self._add_variables(mp)

        # Nodes
        base  = mp.CreateNewNode(3, 0.0, 0.0, 0.0)
        node1 = mp.CreateNewNode(1, 1.0, 0.0, 0.0)
        node2 = mp.CreateNewNode(2, 2.0, 0.0, 0.0)

        self._add_dofs(mp)

        # Constrain transverse motion
        for node in mp.Nodes:
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Z)

        # Fix base in x
        base.Fix(KratosMultiphysics.DISPLACEMENT_X)

        # Properties
        prop1 = mp.GetProperties()[1]
        prop2 = mp.CreateNewProperties(2)

        prop1.SetValue(KratosMultiphysics.YOUNG_MODULUS, stiffness)
        prop1.SetValue(StructuralMechanicsApplication.CROSS_AREA, 1.0)
        prop1.SetValue(KratosMultiphysics.DENSITY, 0.0)
        prop1.SetValue(
            KratosMultiphysics.CONSTITUTIVE_LAW,
            StructuralMechanicsApplication.TrussConstitutiveLaw()
        )

        prop2.SetValue(KratosMultiphysics.YOUNG_MODULUS, stiffness / 2.0)
        prop2.SetValue(StructuralMechanicsApplication.CROSS_AREA, 1.0)
        prop2.SetValue(KratosMultiphysics.DENSITY, 0.0)
        prop2.SetValue(
            KratosMultiphysics.CONSTITUTIVE_LAW,
            StructuralMechanicsApplication.TrussConstitutiveLaw()
        )

        # Elements
        mp.CreateNewElement("LinearTrussElement3D2N", 1, [3, 1], prop1)
        mp.CreateNewElement("LinearTrussElement3D2N", 2, [1, 2], prop2)
        mp.CreateNewElement("NodalConcentratedElement3D1N", 3, [1], prop1)
        mp.CreateNewElement("NodalConcentratedElement3D1N", 4, [2], prop1)

        # Conditions
        cond_real_1 = mp.CreateNewCondition("PointLoadCondition3D1N", 1, [1], prop1)
        cond_real_2 = mp.CreateNewCondition("PointLoadCondition3D1N", 2, [2], prop1)

        cond_imag_1 = mp.CreateNewCondition("ImaginaryPointLoadCondition3D1N", 3, [1], prop1)
        cond_imag_2 = mp.CreateNewCondition("ImaginaryPointLoadCondition3D1N", 4, [2], prop1)

        # Submodelparts for RHS assembly
        real_loads_1_smp = mp.CreateSubModelPart("RealLoads1")
        real_loads_2_smp = mp.CreateSubModelPart("RealLoads2")
        imag_loads_smp = mp.CreateSubModelPart("ImaginaryLoads")

        # Add nodes used by the conditions
        real_loads_1_smp.AddNode(node1, 0)
        real_loads_2_smp.AddNode(node2, 0)
        imag_loads_smp.AddNode(node1, 0)
        imag_loads_smp.AddNode(node2, 0)

        # Add conditions
        real_loads_1_smp.AddCondition(cond_real_1, 0)
        real_loads_2_smp.AddCondition(cond_real_2, 0)
        imag_loads_smp.AddCondition(cond_imag_1, 0)
        imag_loads_smp.AddCondition(cond_imag_2, 0)

        # Lumped masses
        mp.Elements[3].SetValue(KratosMultiphysics.NODAL_MASS, mass)
        mp.Elements[4].SetValue(KratosMultiphysics.NODAL_MASS, mass / 2.0)

        # Real load: node1 only
        node1.SetSolutionStepValue(
            StructuralMechanicsApplication.POINT_LOAD, 0, [1.0, 0.0, 0.0]
        )
        node2.SetSolutionStepValue(
            StructuralMechanicsApplication.POINT_LOAD, 0, [2.0, 0.0, 0.0]
        )

        # Imaginary load: node1 only
        node1.SetSolutionStepValue(
            StructuralMechanicsApplication.POINT_LOAD_IMAGINARY, 0, [1.0, 0.0, 0.0]
        )
        node2.SetSolutionStepValue(
            StructuralMechanicsApplication.POINT_LOAD_IMAGINARY, 0, [0.0, 0.0, 0.0]
        )

        return mp

    def _setup_direct_harmonic_solver(self, mp, echo=0):
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(
            KratosMultiphysics.LinearSolver()
        )

        scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()

        complex_solver_settings = KratosMultiphysics.Parameters(r'''{
            "solver_type" : "skyline_lu_complex"
        }''')

        factory = KratosMultiphysics.ComplexLinearSolverFactory()
        complex_linear_solver = factory.Create(complex_solver_settings)

        direct_harmonic_settings = KratosMultiphysics.Parameters(r'''{
            "assemble_damping_matrix"            : false,
            "real_load_sub_model_part_list"      : [
                "RealLoads1",
                "RealLoads2"
            ],
            "imaginary_load_sub_model_part_list" : [
                "ImaginaryLoads"
            ]
        }''')

        direct_strategy = StructuralMechanicsApplication.DirectHarmonicAnalysisStrategy(
            mp,
            scheme,
            builder_and_solver,
            complex_linear_solver,
            direct_harmonic_settings
        )
        direct_strategy.SetEchoLevel(echo)
        return direct_strategy

    def _compute_reference_solution(self, stiffness, mass, omega):
        k1 = stiffness
        k2 = stiffness / 2.0
        m1 = mass
        m2 = mass / 2.0

        K = np.array([
            [k1 + k2, -k2],
            [-k2,      k2]
        ], dtype=float)

        M = np.array([
            [m1, 0.0],
            [0.0, m2]
        ], dtype=float)

        A = K - (omega ** 2) * M
        f = np.array([1.0 + 1.0j, 2.0 + 0.0j], dtype=complex)

        return np.linalg.solve(A.astype(complex), f)

    def test_direct_mdof_harmonic_undamped(self):
        current_model = KratosMultiphysics.Model()

        stiffness = 10.0
        mass = 2.0

        mp = self._create_3dof_geometry(current_model, stiffness, mass)

        test_frequencies = [1.0, 1.2, 1.4]

        for omega in test_frequencies:
            direct_solver = self._setup_direct_harmonic_solver(mp)
            mp.CloneTimeStep(omega)
            direct_solver.Solve()

            u_ref = self._compute_reference_solution(stiffness, mass, omega)

            self.assertAlmostEqual(
                mp.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0),
                u_ref[0].real,
                delta=1e-8
            )
            self.assertAlmostEqual(
                mp.Nodes[2].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0),
                u_ref[1].real,
                delta=1e-8
            )

            self.assertAlmostEqual(
                mp.Nodes[1].GetSolutionStepValue(
                    StructuralMechanicsApplication.DISPLACEMENT_IMAGINARY
                )[0],
                u_ref[0].imag,
                delta=1e-8
            )
            self.assertAlmostEqual(
                mp.Nodes[2].GetSolutionStepValue(
                    StructuralMechanicsApplication.DISPLACEMENT_IMAGINARY
                )[0],
                u_ref[1].imag,
                delta=1e-8
            )


if __name__ == "__main__":
    KratosUnittest.main()