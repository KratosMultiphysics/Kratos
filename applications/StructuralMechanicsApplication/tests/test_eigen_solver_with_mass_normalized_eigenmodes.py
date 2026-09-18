import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import KratosMultiphysics.KratosUnittest as KratosUnittest

# External imports
import os
import json

'''
Test description:
This test does an eigenvalue analysis on a simple cantilever beam, and compares the mass
normalized eigenmodes with a reference solution
'''

def GetFilePath(file_name):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), file_name)

@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class TestEigenSolverWithMassNormalizedEigenmodes(KratosUnittest.TestCase):
    # muting the output
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

    def test_eigen_with_mass_normalized_eigenmodes(self):
        self.execute_test_eigen_with_mass_normalized_eigenmodes(use_block_builder=False)

    def execute_test_eigen_with_mass_normalized_eigenmodes(self, use_block_builder):
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
                "normalize_eigenvectors_with_mass_matrix": true,
                "rotation_dofs"            : true,
                "block_builder"            : true
            }
        }""")

        analysis_parameters["solver_settings"]["block_builder"].SetBool(use_block_builder)
        
        # Create the model, the analysis and run the eigenvalue problem
        model = KratosMultiphysics.Model()
        analysis = StructuralMechanicsAnalysis(model, analysis_parameters.Clone())
        model_part = model["Structure"]
        SetupSystem(model_part)
        analysis.Run()

        # Get the reference file name
        reference_file = self.__GetReferenceFileName()

        # Compare the solution with the reference solution
        self.__CompareEigenSolutionWithReference(model_part, reference_file)
    
    def __GetReferenceFileName(self):
        file_name = os.path.join("eigen_test", "eigen_with_mass_normalized_eigenmodes_reference.json")
        return GetFilePath(file_name)

    """Compare current eigen solution against the stored reference solution."""
    def __CompareEigenSolutionWithReference(self, model_part, file_name):
        with open(file_name, "r") as f:
            ref_data = json.load(f)

        # Compare the eigenvalues
        eigen_val_vec = model_part.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR]
        ref_eigenvalues = ref_data["eigenvalues"]

        self.assertVectorAlmostEqual(eigen_val_vec, ref_eigenvalues)

        # Compare the eigenvectors
        for node in model_part.Nodes:
            node_id_str = str(node.Id)
            self.assertIn(node_id_str, ref_data["nodes"])

            eig_vec_mat = node[StructuralMechanicsApplication.EIGENVECTOR_MATRIX]
            ref_mat_list = ref_data["nodes"][node_id_str]

            rows = len(ref_mat_list)
            cols = len(ref_mat_list[0]) if rows > 0 else 0

            # Build a Kratos Matrix from the reference solution
            ref_mat = KratosMultiphysics.Matrix(rows, cols)
            for i in range(rows):
                for j in range(cols):
                    ref_mat[i, j] = ref_mat_list[i][j]

            self.assertMatrixAlmostEqual(eig_vec_mat, ref_mat)

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
