import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtils
import numpy as np

class TestVectorizedInterpolation(UnitTest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        self.mp = self.model.CreateModelPart("FluidModelPart")
        self.mp.SetBufferSize(2)
        self.mp.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 3)
        self.mp.AddNodalSolutionStepVariable(Kratos.DISTANCE)
        self.mp.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        Kratos.ModelPartIO("./TwoFluidMassConservationProcTest/3Dtest").ReadModelPart(self.mp)

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

        nnodes = len(self.mp.Nodes)
        vdata = np.array(Kratos.VariableUtils().GetSolutionStepValuesVector(self.mp.Nodes, Kratos.VELOCITY, 0, 3)).reshape(nnodes,3)
        coords = np.array(Kratos.VariableUtils().GetCurrentPositionsVector(self.mp.Nodes,3)).reshape(nnodes,3)

        xyz = np.array([
            [0.2,0.3,0.4],
            [0.1,0.1,0.3]
        ])
        print("xyz",xyz)

        elem_ids, geom_ids,N = locator.VectorizedFind(xyz)
        N = np.array(N)
        geom_ids = np.array(geom_ids, dtype=int)
        geom_ids -= 1 #indices are IDs-1


        interpolated_data = self.Interpolate(geom_ids, N, vdata)
        print("interpolated",interpolated_data)



    def tearDown(self):
        KratosUtils.DeleteFileIfExisting("Cavity/square5.time")


if __name__ == '__main__':
    UnitTest.main()
