import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import KratosMultiphysics.KratosUnittest as KratosUnittest

'''
Test description:
This test does an eigenvalue analysis on a simple cantilever beam
It runs the test in two configurations:
1. Regular setup of the system
2. One node in the middle is duplicated. The two sides of the cantilever are connected with constraints
The test checks if the two solutions are the same, which should be the case if everything works fine
'''

class StructuralMechanicsAnalysisWithConstraints(StructuralMechanicsAnalysis):
    def ModifyInitialGeometry(self):
        super(StructuralMechanicsAnalysisWithConstraints, self).ModifyInitialGeometry()

        model_part = self.model["Structure"]
        num_nodes = model_part.NumberOfNodes()

        constraint_node_id = int(num_nodes/2)
        constraint_node = model_part.Nodes[constraint_node_id]
        aux_node_id = num_nodes # note that this is different from before bcs now there is also the constraint node in the model-part

        aux_node = model_part.Nodes[aux_node_id]

        vars_to_constrain = [
            KratosMultiphysics.DISPLACEMENT_X,
            KratosMultiphysics.DISPLACEMENT_Y,
            KratosMultiphysics.DISPLACEMENT_Z,
            KratosMultiphysics.ROTATION_X,
            KratosMultiphysics.ROTATION_Y,
            KratosMultiphysics.ROTATION_Z
        ]

        for i_var, var in enumerate(vars_to_constrain):
            model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", i_var+1, constraint_node, var, aux_node, var, 1.0, 0.0)


@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class TestEigenSolverWithConstraints(KratosUnittest.TestCase):
    # muting the output
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

    def test_eigen_with_constraints_block_builder(self):
        self.execute_test_eigen_with_constraints(use_block_builder=True)

    def test_eigen_with_constraints_elimination_builder(self):
        self.execute_test_eigen_with_constraints(use_block_builder=False)

    def execute_test_eigen_with_constraints(self, use_block_builder):
        analysis_parameters = KratosMultiphysics.Parameters("""{
            "problem_data"    : {
                "parallel_type" : "OpenMP",
                "echo_level"    : 1,
                "start_time"    : 0.0,
                "end_time"      : 1.0
            },
            "solver_settings" : {
                "solver_type"              : "eigen_value",
                "model_part_name"          : "Structure",
                "domain_size"              : 3,
                "model_import_settings"    : {
                    "input_type"     : "use_input_model_part"
                },
                "time_stepping"            : {
                    "time_step" : 1.1
                },
                "use_computing_model_part" : false,
                "rotation_dofs"            : true,
                "block_builder"            : true
            }
        }""")

        analysis_parameters["solver_settings"]["block_builder"].SetBool(use_block_builder)
        analysis_parameters_with_constraints = analysis_parameters.Clone()
        analysis_parameters_with_constraints["solver_settings"]["block_builder"].SetBool(True) # Currently the EliminationB&S does not reliably work with constraints

        model = KratosMultiphysics.Model()
        analysis = StructuralMechanicsAnalysis(model, analysis_parameters.Clone())
        model_part = model["Structure"]
        SetupSystem(model_part, use_constraints=False)
        analysis.Run()

        model_with_constraints = KratosMultiphysics.Model()
        analysis_with_constraints = StructuralMechanicsAnalysisWithConstraints(model_with_constraints, analysis_parameters_with_constraints)
        model_part_with_constraints = model_with_constraints["Structure"]
        SetupSystem(model_part_with_constraints, use_constraints=True)
        analysis_with_constraints.Run()

        self.__CompareEigenSolution(model_part, model_part_with_constraints)

    def __CompareEigenSolution(self, model_part, model_part_with_constraints):
        eigen_val_vec = model_part.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR]
        eigen_val_vec_with_constraints = model_part_with_constraints.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR]

        self.assertEqual(eigen_val_vec.Size(), eigen_val_vec_with_constraints.Size())

        for eig_val, eig_val_constr in zip(eigen_val_vec, eigen_val_vec_with_constraints):
            self.assertAlmostEqual(eig_val, eig_val_constr, 8)

        for node in model_part.Nodes:
            node_const = model_part_with_constraints.Nodes[node.Id] # to make sure to get the corresponding node
            eig_vec_mat = node[StructuralMechanicsApplication.EIGENVECTOR_MATRIX]
            eig_vec_mat_contr = node_const[StructuralMechanicsApplication.EIGENVECTOR_MATRIX]

            self.__CompareMatrix(eig_vec_mat, eig_vec_mat_contr, 10) # Note: this might me too strict depending on the eigenvalue solver (works fine with eigen_eigensystem in compination with the eigen sparse-lu)

    def __CompareMatrix(self, mat_1, mat_2, tol=7):
        self.assertEqual(mat_1.Size1(), mat_2.Size1())
        self.assertEqual(mat_1.Size2(), mat_2.Size2())

        for i in range(mat_1.Size1()):
            for j in range(mat_1.Size2()):
                self.assertAlmostEqual(abs(mat_1[i,j]), abs(mat_2[i,j]), tol)


def SetupSystem(model_part, use_constraints):
    # manually creating the system to avoid extra files and because this way the manipulation of the system for the constraints is easier
    num_nodes = 10
    for i in range(num_nodes):
        model_part.CreateNewNode(i+1, 1.0*i, 0.0, 0.0)

    clamped_node = model_part.Nodes[1]

    # apply cantilever boundary conditions
    clamped_node.Fix(KratosMultiphysics.DISPLACEMENT_X)
    clamped_node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
    clamped_node.Fix(KratosMultiphysics.DISPLACEMENT_Z)

    clamped_node.Fix(KratosMultiphysics.ROTATION_X)
    clamped_node.Fix(KratosMultiphysics.ROTATION_Y)
    clamped_node.Fix(KratosMultiphysics.ROTATION_Z)

    element_connectivities = []
    for i in range(num_nodes-1):
        element_connectivities.append([i+1, i+2])

    if use_constraints:
        # create a duplicated node to apply the constraints
        constraint_node_id = int(num_nodes/2)
        constraint_node = model_part.Nodes[constraint_node_id]
        aux_node_id = num_nodes+1
        model_part.CreateNewNode(aux_node_id, constraint_node.X, 0.0, 0.0)

        # modify the connectivities to use the aux node
        element_connectivities[constraint_node_id-1][0] = aux_node_id

    props = model_part.CreateNewProperties(0)
    props[KratosMultiphysics.YOUNG_MODULUS] = 210e9
    props[KratosMultiphysics.DENSITY] = 7850
    props[StructuralMechanicsApplication.CROSS_AREA] = 0.01
    props[KratosMultiphysics.POISSON_RATIO] = 0.30
    props[StructuralMechanicsApplication.TORSIONAL_INERTIA] = 0.00001
    props[StructuralMechanicsApplication.I22] = 0.00002
    props[StructuralMechanicsApplication.I33] = 0.00001

    for i_elem, connectivity in enumerate(element_connectivities):
        model_part.CreateNewElement("CrLinearBeamElement3D2N", i_elem+1, connectivity, props)


if __name__ == '__main__':
    KratosUnittest.main()
