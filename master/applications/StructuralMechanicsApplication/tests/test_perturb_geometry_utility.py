
import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.StructuralMechanicsApplication.perturb_geometry_sparse_utility import PerturbGeometrySparseUtility
from KratosMultiphysics.StructuralMechanicsApplication.perturb_geometry_subgrid_utility import PerturbGeometrySubgridUtility

import math

class SparseUtilityCustom(PerturbGeometrySparseUtility):
    """ SparseUtilityCustom
    This class is derived to override the PerturbGeometry method
    """
    def PerturbGeometry(self, mp):
        # Apply random field vectors to geometry
        self.utility.ApplyRandomFieldVectorsToGeometry(mp, [1,0,0,0,0])

class SubgridUtilityCustom(PerturbGeometrySubgridUtility):
    """SubgridUtilityCustom
    This class is derived to override the PerturbGeometry method
    """
    def PerturbGeometry(self, mp):
        # Apply random field vectors to geometry
        self.utility.ApplyRandomFieldVectorsToGeometry(mp, [1,0,0,0,0])

# This test generates a random field for a square plate with 5x5 nodes with the sparse and the subgrid method.
# The test is passed when the first perturbation vectors from both models are equal.
class BaseTestPerturbGeometryUtility(KratosUnittest.TestCase):
    """BaseTestPerturbGeometryUtility
    Base class of the test
    """
    @classmethod
    def _add_dofs(cls, mp):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z, mp)

    @classmethod
    def _add_variables(cls, mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

    @classmethod
    def _create_nodes(cls, mp, NumOfNodes, length):
        counter = 0
        for y in range(NumOfNodes):
            for x in range(NumOfNodes):
                counter = counter + 1
                mp.CreateNewNode(counter, x/(NumOfNodes-1)*length, y/(NumOfNodes-1)*length,0.0)

    @classmethod
    def _create_elements(cls, mp, NumOfNodes):
        element_name = "ShellThinElementCorotational3D4N"
        counter = 0
        for y in range( NumOfNodes-1 ):
            for x in range(NumOfNodes - 1):
                # Aligned counter-clockwise
                counter = counter + 1
                node1 = NumOfNodes*(y)+(x+1)
                node2 = NumOfNodes*(y)+(x+2)
                node3 = NumOfNodes*(y+1) + (x+2)
                node4 = NumOfNodes*(y+1) + (x+1)
                mp.CreateNewElement( element_name, counter, [ node1, node2, node3, node4 ],mp.GetProperties()[0] )

    def _set_up_system(self, model, NumOfNodes, length):
        mp = model.CreateModelPart("Structure")
        self._add_variables(mp)
        self._create_nodes(mp,NumOfNodes, length)
        self._add_dofs(mp)
        self._create_elements(mp, NumOfNodes)
        return mp

    def _compare_random_field_vectors(self, mp1, mp2):
        nodes1 = mp1.GetNodes()
        nodes2 = mp2.GetNodes()
        sum_ = 0
        for node1, node2 in zip(nodes1, nodes2):
            sum_ += (node1.Z0 - node2.Z0)**2

        self.assertLess( math.sqrt(sum_), 1.0e-10)

@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class TestPerturbGeometryUtility(BaseTestPerturbGeometryUtility):
    def test_perturb_geometry_utility(self):
        num_of_nodes_per_egde = 5
        length = 1000
        # Sparse method
        model_sparse =  KratosMultiphysics.Model()
        mp_sparse = self._set_up_system(model_sparse, num_of_nodes_per_egde, length)
        settings = KratosMultiphysics.Parameters("""{
                    "eigensolver_settings"  : {
                        "solver_type"               : "eigen_eigensystem",
                        "max_iteration"             : 1000,
                        "tolerance"                 : 1e-10,
                        "number_of_eigenvalues"     : 10,
                        "normalize_eigenvectors"    : false,
                        "echo_level"                : 0,
                        "use_mkl_if_available"      : false
                        },
                    "perturbation_settings" : {
                        "max_displacement"          : 1,
                        "correlation_length"        : 500,
                        "truncation_error"          : 0.1,
                        "echo_level"                : 0
                    }
                }""")
        SparseUtility = SparseUtilityCustom(mp_sparse, settings)
        SparseUtility.PerturbGeometry(mp_sparse)
        # Subgrid method
        model_subgrid =  KratosMultiphysics.Model()
        mp_subgrid = self._set_up_system(model_subgrid, num_of_nodes_per_egde, length)
        settings = KratosMultiphysics.Parameters("""{
            "eigensolver_settings"  : {
                "solver_type"               : "dense_eigensolver",
                "ascending_order"           : false,
                "echo_level"                : 0
                },
            "perturbation_settings" : {
                "min_distance_subgrid"      : 1,
                "max_displacement"          : 1.0,
                "correlation_length"        : 500,
                "truncation_error"          : 0.2,
                "echo_level"                : 0
            }
        }""")
        SubgridUtility = SubgridUtilityCustom(mp_subgrid, settings)
        SubgridUtility.PerturbGeometry(mp_subgrid)

        # Check if first random field vectors are equal
        self._compare_random_field_vectors(mp_sparse, mp_subgrid)

if __name__ == '__main__':
    KratosUnittest.main()