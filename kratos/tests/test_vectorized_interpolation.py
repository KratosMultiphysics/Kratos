import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtils
import numpy as np
import os

class TestVectorizedInterpolation(UnitTest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        self.mp = self.model.CreateModelPart("FluidModelPart")
        self.mp.SetBufferSize(2)
        self.mp.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 3)
        self.mp.AddNodalSolutionStepVariable(Kratos.DISTANCE)
        self.mp.AddNodalSolutionStepVariable(Kratos.VELOCITY)

        Kratos.ModelPartIO(os.path.join(".","auxiliar_files_for_python_unittest","mdpa_files","coarse_sphere")).ReadModelPart(self.mp)

        #initialize VELOCITY on nodes with the spatial position. This will allow checking analitically for the interpolation
        for node in self.mp.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY,0,[node.X,node.Y,node.Z])

    def Interpolate(self,node_indices, N, origin_data):
        data_gathered = origin_data[node_indices]

        # 2. Expand dimensions of N to allow broadcasting for the multiplication
        #    Resulting shape will be (number_of_particles, number_of_nodes_per_particle, 1)
        N_expanded = N[..., np.newaxis]

        # 3. Multiply and sum over the second axis (number_of_nodes_per_particle)
        #    Final shape will be (number_of_particles, domain_size)
        vp = np.sum(N_expanded * data_gathered, axis=1)

        return vp

    def testVectorizedInterpolation(self):
        locator = Kratos.BinBasedFastPointLocator3D(self.mp)
        locator.UpdateSearchDatabase()

        #obtaining a copy of database onto numpy array, and reshaping it into matrices
        nnodes = len(self.mp.Nodes)
        vdata = np.array(Kratos.VariableUtils().GetSolutionStepValuesVector(self.mp.Nodes, Kratos.VELOCITY, 0, 3),copy=False).reshape(nnodes,3)
        coords = np.array(Kratos.VariableUtils().GetCurrentPositionsVector(self.mp.Nodes,3)).reshape(nnodes,3)

        ##coordinates to search for
        xyz = np.array([
            [0.02,0.03,0.04],
            [0.07,0.01,0.03]
        ])
        print("50")
        elem_ids, geom_ids,N = locator.VectorizedFind(xyz)
        print(elem_ids)
        #### here we are transofrming Ids into indices to address vdata and coords.
        #### In the current example it simply involves adding -1 as the Kratos Ids start on 1
        #### on more complex examples, where the nodes in modelpart are not numbered consecutively, this would be a more complex op
        geom_ids -= 1

        #### calling pure-python function for interpolation (vectorized)
        interpolated_data = self.Interpolate(geom_ids, N, vdata)

        ##check for the result of interpolation
        err_norm = np.linalg.norm(xyz - interpolated_data)
        self.assertAlmostEqual(err_norm, 0.0)

if __name__ == '__main__':
    UnitTest.main()
