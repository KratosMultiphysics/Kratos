import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import KratosMultiphysics.KratosUnittest as KratosUnittest

try:
    import scipy
    scipy_is_available = True
except ImportError:
    scipy_is_available = False
'''
Test description:
This test does an eigenvalue analysis on a simple cantilever beam
It compares the results between the scipy solver and the standard eigen solver.
'''

@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
@KratosUnittest.skipUnless(scipy_is_available,"Missing required module: Scipy")
class TestCustomScipyBaseSolver(KratosUnittest.TestCase):
    # muting the output
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

    def test_custom_scipy_base_solver(self):
        analysis_parameters_scipy = KratosMultiphysics.Parameters("""{
            "problem_data"    : {
                "parallel_type" : "OpenMP",
                "echo_level"    : 1,
                "start_time"    : 0.0,
                "end_time"      : 1.0
            },
            "solver_settings" : {
                "solver_type"              : "KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_custom_scipy_base_solver",
                "model_part_name"          : "Structure",
                "domain_size"              : 3,
                "model_import_settings"    : {
                    "input_type"     : "use_input_model_part"
                },
                "time_stepping"            : {
                    "time_step" : 1.1
                },
                "rotation_dofs"            : true,
                "builder_and_solver_settings" : {
                    "use_block_builder" : false
                }
            }
        }""")

        analysis_parameters_eigen = KratosMultiphysics.Parameters("""{
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
                "rotation_dofs"            : true,
                "eigensolver_settings"     : {
                   "solver_type"           : "eigen_eigensystem",
                   "number_of_eigenvalues" : 5,
                   "normalize_eigenvectors": true,
                   "max_iteration"         : 10000,
                   "tolerance"             : 1e-6
                },
                "builder_and_solver_settings" : {
                    "use_block_builder" : true
                }
            }
        }""")

        model_scipy = KratosMultiphysics.Model()
        analysis_scipy = StructuralMechanicsAnalysis(model_scipy, analysis_parameters_scipy.Clone())
        model_part_scipy = model_scipy["Structure"]
        SetupSystem(model_part_scipy)
        analysis_scipy.Run()

        model_eigen = KratosMultiphysics.Model()
        analysis_eigen = StructuralMechanicsAnalysis(model_eigen, analysis_parameters_eigen.Clone())
        model_part_eigen = model_eigen["Structure"]
        SetupSystem(model_part_eigen)
        analysis_eigen.Run()

        self.__CompareEigenSolution(model_part_scipy, model_part_eigen)

    def __CompareEigenSolution(self, model_part, model_part_ref):
        eigen_val_vec = model_part.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR]
        eigen_val_vec_ref = model_part_ref.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR]

        self.assertEqual(eigen_val_vec.Size(), eigen_val_vec_ref.Size())

        for eig_val, eig_val_constr in zip(eigen_val_vec, eigen_val_vec_ref):
            self.assertAlmostEqual(eig_val, eig_val_constr, 6)

        for node in model_part.Nodes:
            node_ref = model_part_ref.Nodes[node.Id] # to make sure to get the corresponding node
            eig_vec_mat = node[StructuralMechanicsApplication.EIGENVECTOR_MATRIX]
            eig_vec_mat_ref = node_ref[StructuralMechanicsApplication.EIGENVECTOR_MATRIX]
            self.__CompareMatrix(eig_vec_mat, eig_vec_mat_ref, 8)

    def __CompareMatrix(self, mat_1, mat_2, tol=7):
        self.assertEqual(mat_1.Size1(), mat_2.Size1())
        self.assertEqual(mat_1.Size2(), mat_2.Size2())
        for i in range(mat_1.Size1()):
            for j in range(mat_1.Size2()):
                self.assertAlmostEqual(abs(mat_1[i,j]), abs(mat_2[i,j]), tol)


def SetupSystem(model_part):
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

    props = model_part.CreateNewProperties(0)
    props[KratosMultiphysics.YOUNG_MODULUS] = 210e9
    props[KratosMultiphysics.DENSITY] = 7850
    props[StructuralMechanicsApplication.CROSS_AREA] = 0.01
    props[KratosMultiphysics.POISSON_RATIO] = 0.30
    props[StructuralMechanicsApplication.TORSIONAL_INERTIA] = 0.00001
    props[StructuralMechanicsApplication.I22] = 0.00002
    props[StructuralMechanicsApplication.I33] = 0.00001
    props[KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX] = False

    for i_elem, connectivity in enumerate(element_connectivities):
        model_part.CreateNewElement("CrLinearBeamElement3D2N", i_elem+1, connectivity, props)


if __name__ == '__main__':
    KratosUnittest.main()
