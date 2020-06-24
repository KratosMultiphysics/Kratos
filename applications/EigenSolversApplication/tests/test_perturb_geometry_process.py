from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

try:
    import KratosMultiphysics.EigenSolversApplication as EigenSolversApplication
    eigen_solvers_is_available = True
except ImportError:
    eigen_solvers_is_available = False

import numpy as np

# This test generates a random field for a sqaure plate with 5x5 nodes with the sparse and the subgrid method
# The test is passed when the first perturbation vectors from both models are equal
class BaseTestPerturbGeometryProcess(KratosUnittest.TestCase):
    @classmethod
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
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

    def _create_nodes(self,mp, NumOfNodes,length):
        # Create nodes
        counter = 0
        for y in range(NumOfNodes):
            for x in range(NumOfNodes):
                counter = counter + 1
                mp.CreateNewNode(counter,x/(NumOfNodes-1)*length,y/(NumOfNodes-1)*length,0.0)

    def _create_elements(self,mp,NumOfNodes):
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

    def _set_up_system(self,model,NumOfNodes,length):
        mp = model.CreateModelPart("Structure")
        self._add_variables(mp)
        self._create_nodes(mp,NumOfNodes,length)
        self._add_dofs(mp)
        self._create_elements(mp,NumOfNodes)
        return mp

    def _compare_random_field_vectors(self, mp1, mp2):
        nodes1 = mp1.GetNodes()
        nodes2 = mp2.GetNodes()
        sum = 0
        for node1, node2 in zip(nodes1,nodes2):
            sum += (node1.Z0 - node2.Z0)**2

        self.assertLess( np.sqrt(sum), 1.0e-10)

class TestPerturbGeometryProcess(BaseTestPerturbGeometryProcess):
    @KratosUnittest.skipUnless(eigen_solvers_is_available,"EigenSolversApplication not available")
    def test_perturb_geometry_process(self):
        num_of_nodes_per_egde = 5
        length = 1000
        # Subgrid method
        model_subgrid =  KratosMultiphysics.Model()
        mp_subgrid = self._set_up_system(model_subgrid ,num_of_nodes_per_egde,length )
        subgrid_method = EigenSolversApplication.PerturbGeometrySubgridProcess(mp_subgrid, 1.0)
        subgrid_method.SetEchoLevel(0)
        minDistance = 1
        correlation_length = 500;
        truncation_tolerance = 0.2
        num_variables_subgrid = subgrid_method.CreateEigenvectors( mp_subgrid, minDistance, correlation_length, truncation_tolerance)
        subgrid_method.AssembleEigenvectors( mp_subgrid, [0,0,0,0,1] )
        # Sparse method
        model_sparse =  KratosMultiphysics.Model()
        mp_sparse = self._set_up_system(model_sparse ,num_of_nodes_per_egde,length )
        eigensolver_settings = KratosMultiphysics.Parameters("""
            {
                "max_iteration"         : 1000,
                "tolerance"             : 1e-6,
                "number_of_eigenvalues" : 10,
                "echo_level"            : 0,
                "normalize_eigenvectors": false
            }
            """)
        sparse_method = KratosMultiphysics.EigenSolversApplication.PerturbGeometrySparseProcess(mp_sparse,1.0)
        sparse_method.SetEchoLevel(0)
        correlation_length = 500;
        truncation_tolerance = 0.1
        num_variables_sparse = sparse_method.CreateEigenvectors( mp_sparse, correlation_length, truncation_tolerance, eigensolver_settings)
        sparse_method.AssembleEigenvectors( mp_sparse, [1,0,0,0,0] )
        # Check if first random field vectors are equal
        self._compare_random_field_vectors(mp_sparse,mp_subgrid)

if __name__ == '__main__':
    KratosUnittest.main()
