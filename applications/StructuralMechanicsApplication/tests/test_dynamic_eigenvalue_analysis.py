import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import kratos_utilities


if kratos_utilities.CheckIfApplicationsAvailable("LinearSolversApplication"):
    from KratosMultiphysics import LinearSolversApplication

#Test of eigenvalue analysis and modal decomposition with beam elements according to:
#C. Petersen, Dynamic der Baukonstruktionen, Viehweg Verlag, 2000, p. 252

class BaseTestDynamicEigenvalueAnalysis(KratosUnittest.TestCase):
    def _add_dofs(self,mp):
        # Adding dofs AND their corresponding reactions
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,mp)


    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)


    def _create_nodes(self,mp):
        mp.CreateNewNode(1,0.0000000000,12.0000000000,0.0000000000)
        mp.CreateNewNode(2,0.0000000000,10.8000000000,0.0000000000)
        mp.CreateNewNode(3,0.0000000000,9.6000000000,0.0000000000)
        mp.CreateNewNode(4,0.0000000000,8.4000000000,0.0000000000)
        mp.CreateNewNode(5,0.0000000000,7.2000000000,0.0000000000)
        mp.CreateNewNode(6,0.0000000000,6.0000000000,0.0000000000)
        mp.CreateNewNode(7,0.0000000000,4.8000000000,0.0000000000)
        mp.CreateNewNode(8,0.0000000000,3.6000000000,0.0000000000)
        mp.CreateNewNode(9,0.0000000000,2.4000000000,0.0000000000)
        mp.CreateNewNode(10,0.0000000000,1.2000000000,0.0000000000)
        mp.CreateNewNode(11,0.0000000000,0.0000000000,0.0000000000)

    def _create_elements(self,mp):
        element_name = "CrLinearBeamElement2D2N"
        mp.CreateNewElement(element_name,1,[11,10],mp.GetProperties()[0])
        mp.CreateNewElement(element_name,2,[10, 9],mp.GetProperties()[0])
        mp.CreateNewElement(element_name,3,[9, 8],mp.GetProperties()[0])
        mp.CreateNewElement(element_name,4,[8, 7],mp.GetProperties()[0])
        mp.CreateNewElement(element_name,5,[7, 6],mp.GetProperties()[0])
        mp.CreateNewElement(element_name,6,[6, 5],mp.GetProperties()[0])
        mp.CreateNewElement(element_name,7,[5, 4],mp.GetProperties()[0])
        mp.CreateNewElement(element_name,8,[4, 3],mp.GetProperties()[0])
        mp.CreateNewElement(element_name,9,[3, 2],mp.GetProperties()[0])
        mp.CreateNewElement(element_name,10,[2, 1],mp.GetProperties()[0])

        mass_1 = mp.CreateNewElement("NodalConcentratedElement2D1N", 11, [1], mp.GetProperties()[0])
        mass_2 = mp.CreateNewElement("NodalConcentratedElement2D1N", 12, [6], mp.GetProperties()[0])
        mass_1.SetValue(KratosMultiphysics.NODAL_MASS,250)
        mass_2.SetValue(KratosMultiphysics.NODAL_MASS,500)
        mass_1.Initialize(mp.ProcessInfo)
        mass_2.Initialize(mp.ProcessInfo)

    def _apply_material_properties(self,mp):
        #define properties
        mp.GetProperties()[0].SetValue(KratosMultiphysics.YOUNG_MODULUS,2e7)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.DENSITY,0.01)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.CROSS_AREA,10)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.I33,1.0)

        #g = [0,0,0]
        #mp.GetProperties()[0].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)

        cl = StructuralMechanicsApplication.LinearElastic3DLaw()
        mp.GetProperties()[0].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

    def _apply_BCs(self,mp):
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.ROTATION_Z, True, mp.Nodes)

    def _solve_eigenvalue_problem(self, mp, use_block_builder, echo=0):
        eigensolver_settings = KratosMultiphysics.Parameters("""
        {
            "max_iteration"         : 1000,
            "tolerance"             : 1e-6,
            "number_of_eigenvalues" : 2,
            "echo_level"            : 0,
            "normalize_eigenvectors": true
        }
        """)

        eigen_solver = LinearSolversApplication.EigensystemSolver(eigensolver_settings)
        if use_block_builder:
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(eigen_solver)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(eigen_solver)

        eigen_scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
        compute_modal_decomposition = True
        # see "structural_mechanics_eigensolver.py", values are for the "EigensystemSolver"
        mass_matrix_diagonal_value = 0.0
        stiffness_matrix_diagonal_value = 1.0
        eig_strategy = StructuralMechanicsApplication.EigensolverStrategy(mp,
                                                                          eigen_scheme,
                                                                          builder_and_solver,
                                                                          mass_matrix_diagonal_value,
                                                                          stiffness_matrix_diagonal_value,
                                                                          compute_modal_decomposition)
        eig_strategy.SetEchoLevel(echo)
        eig_strategy.Solve()


    def _set_up_system(self,current_model):
        mp = current_model.CreateModelPart("Structure")
        #mp.SetBufferSize(2)

        self._add_variables(mp)
        self._apply_material_properties(mp)
        self._create_nodes(mp)
        self._add_dofs(mp)
        self._create_elements(mp)

        #create a submodelpart for dirichlet boundary conditions
        bcs_dirichlet = mp.CreateSubModelPart("BoundaryCondtionsDirichlet")
        bcs_dirichlet.AddNodes([11])
        self._apply_BCs(bcs_dirichlet)

        return mp

    def _check_eigenvalue_results(self,mp,reference_eigenvalues):
        #check that the results are exact on the node
        eigenvalues = mp.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR]
        self.assertAlmostEqual(eigenvalues[0], reference_eigenvalues[0], 4)
        self.assertAlmostEqual(eigenvalues[1], reference_eigenvalues[1], 4)

    def _check_modal_decomposition_results(self,mp):
        #check that the results are exact on the node
        eigenvalues = mp.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR]
        modal_mass = mp.ProcessInfo[StructuralMechanicsApplication.MODAL_MASS_MATRIX]
        modal_stiffness = mp.ProcessInfo[StructuralMechanicsApplication.MODAL_STIFFNESS_MATRIX]
        self.assertAlmostEqual(modal_mass[0,0], 1.0, 4)
        self.assertAlmostEqual(modal_mass[1,1], 1.0, 4)
        self.assertAlmostEqual(modal_mass[1,0], 0.0, 4)
        self.assertAlmostEqual(modal_mass[0,1], 0.0, 4)
        self.assertAlmostEqual(modal_stiffness[0,0], eigenvalues[0], 4)
        self.assertAlmostEqual(modal_stiffness[1,1], eigenvalues[1], 4)
        self.assertAlmostEqual(modal_stiffness[1,0], 0.0, 4)
        self.assertAlmostEqual(modal_stiffness[0,1], 0.0, 4)


@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class TestDynamicEigenvalueAnalysis(BaseTestDynamicEigenvalueAnalysis):

    def test_dynamic_eigenvalue_analysis_block_builder(self):
        self.execute_test_dynamic_eigenvalue_analysis(use_block_builder=True)

    def test_dynamic_eigenvalue_analysis_elimination_builder(self):
        self.execute_test_dynamic_eigenvalue_analysis(use_block_builder=False)

    def execute_test_dynamic_eigenvalue_analysis(self, use_block_builder):
        reference_eigenvalues = [115.1882,3056.9526]
        current_model = KratosMultiphysics.Model()
        mp = self._set_up_system(current_model)
        self._solve_eigenvalue_problem(mp, use_block_builder)

        self._check_eigenvalue_results(mp,reference_eigenvalues)

        self._check_modal_decomposition_results(mp)

if __name__ == '__main__':
    KratosUnittest.main()
