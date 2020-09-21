import KratosMultiphysics
from KratosMultiphysics import KratosUnittest
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
data_comm = KratosMultiphysics.DataCommunicator.GetDefault()

class TestTransferConditions(KratosUnittest.TestCase):

    def _create_nodes(self, current_mp, dimension):
        if dimension is 2:
            current_mp.CreateNewNode(1,0.0,0.0,0.0)
            current_mp.CreateNewNode(2,0.0,1.0,0.0)
            current_mp.CreateNewNode(3,1.0,0.0,0.0)

    def _create_condition(self, mp, dimension, condition_type, id):
        #Create nodes
        self._create_nodes(mp, dimension)

        # Ensure that the property 1 is created
        mp.GetProperties()[1]

        if dimension == 2:
            if condition_type == "dirichlet":
                cond = mp.CreateNewCondition("MPMParticlePenaltyDirichletCondition2D3N", id, [1, 2, 3], mp.GetProperties()[1])
            if condition_type == "neumann":
                cond = mp.CreateNewCondition("MPMParticlePointLoadCondition2D3N", id, [1, 2, 3], mp.GetProperties()[1])

    def _assign_pseudo_variables(self, cond, condition_type):
        SomeProcessInfo = KratosMultiphysics.ProcessInfo()
        normal_vector = [KratosMultiphysics.Vector([1.0,0.0,0.0])]
        cond.SetValuesOnIntegrationPoints(KratosParticle.MPC_NORMAL,normal_vector,SomeProcessInfo)
        xg = [KratosMultiphysics.Vector([1.5,-1.0,2.1])]
        cond.SetValuesOnIntegrationPoints(KratosParticle.MPC_COORD,xg,SomeProcessInfo)
        velocity = [KratosMultiphysics.Vector([1.5,-1.0,2.45])]
        cond.SetValuesOnIntegrationPoints(KratosParticle.MPC_VELOCITY,velocity,SomeProcessInfo)
        acceleration = [KratosMultiphysics.Vector([1.5,-1.12,2.45])]
        cond.SetValuesOnIntegrationPoints(KratosParticle.MPC_ACCELERATION,acceleration,SomeProcessInfo)
        if condition_type == "dirichlet":
            cond.SetValuesOnIntegrationPoints(KratosParticle.PENALTY_FACTOR,[100.0],SomeProcessInfo)
            displacement = [KratosMultiphysics.Vector([1.22,-1.11,0.0])]
            cond.SetValuesOnIntegrationPoints(KratosParticle.MPC_DISPLACEMENT,displacement,SomeProcessInfo)
            imposed_displacement = [KratosMultiphysics.Vector([1.0,-1.0,0.0])]
            cond.SetValuesOnIntegrationPoints(KratosParticle.MPC_IMPOSED_DISPLACEMENT,imposed_displacement,SomeProcessInfo)
            imposed_velocity = [KratosMultiphysics.Vector([1.0,-1.0,1.1])]
            cond.SetValuesOnIntegrationPoints(KratosParticle.MPC_IMPOSED_VELOCITY,imposed_velocity,SomeProcessInfo)
            imposed_acceleration = [KratosMultiphysics.Vector([1.0,-1.0,2.1])]
            cond.SetValuesOnIntegrationPoints(KratosParticle.MPC_IMPOSED_ACCELERATION,imposed_acceleration,SomeProcessInfo)
        else:
            point_load = [KratosMultiphysics.Vector([3.3,4.4,5.5])]
            cond.SetValuesOnIntegrationPoints(KratosParticle.POINT_LOAD,point_load ,SomeProcessInfo)

    def _check_conditions(self, mp, dimension):
        SomeProcessInfo = KratosMultiphysics.ProcessInfo()
        for cond in mp.Conditions:
            #Check geometry
            self.assertEqual( cond.GetGeometry().WorkingSpaceDimension(), dimension)
            self.assertEqual( cond.GetGeometry().LocalSpaceDimension(), dimension)
            self.assertEqual( cond.GetGeometry().IntegrationPointsNumber(), 1)
            jacobian = cond.GetGeometry().Jacobian(0)
            shape_functions_values = cond.GetGeometry().ShapeFunctionsValues()
            shape_functions_derivatives = cond.GetGeometry().ShapeFunctionDerivatives(1,0)
            if dimension is 2:
                jacobian_ref = KratosMultiphysics.Matrix(2,2,0.0)
                jacobian_ref[0,1] = 1.0
                jacobian_ref[1,0] = 1.0
                shape_functions_values_ref = KratosMultiphysics.Matrix(1,3,1.0/3.0)
                shape_functions_derivatives_ref = KratosMultiphysics.Matrix(3,2,0.0)
                shape_functions_derivatives_ref[0,0] = -1.0
                shape_functions_derivatives_ref[1,0] = 1.0
                shape_functions_derivatives_ref[0,1] = -1.0
                shape_functions_derivatives_ref[2,1] = 1.0

                self.assertMatrixAlmostEqual(jacobian, jacobian_ref, 7)
                self.assertMatrixAlmostEqual(shape_functions_values, shape_functions_values_ref, 7)
                self.assertMatrixAlmostEqual(shape_functions_derivatives, shape_functions_derivatives_ref, 7)
                self.assertVectorAlmostEqual( cond.GetGeometry().Center(), [1.0/3.0, 1.0/3.0, 0.0])
            ##Check condition properties
            unit_normal = cond.CalculateOnIntegrationPoints(KratosParticle.MPC_NORMAL, SomeProcessInfo)
            self.assertVectorAlmostEqual(unit_normal[0],[1.0,0.0,0.0],7)
            if(cond.Info() == "Condition #3"):
                #point_load_condition members
                point_load = cond.CalculateOnIntegrationPoints(KratosParticle.POINT_LOAD, SomeProcessInfo)
                self.assertVectorAlmostEqual(point_load[0],[3.3,4.4,5.5])
            else:
                #penalty_dirichlet_condition members
                penalty_factor = cond.CalculateOnIntegrationPoints(KratosParticle.PENALTY_FACTOR, SomeProcessInfo)
                self.assertAlmostEqual(penalty_factor[0], 100.0, 7)
                #base_dirichlet_condition members
                displacement = cond.CalculateOnIntegrationPoints(KratosParticle.MPC_DISPLACEMENT, SomeProcessInfo)
                self.assertVectorAlmostEqual(displacement[0],[1.22,-1.11,0.0],7)
                imposed_displacement = cond.CalculateOnIntegrationPoints(KratosParticle.MPC_IMPOSED_DISPLACEMENT, SomeProcessInfo)
                self.assertVectorAlmostEqual(imposed_displacement[0],[1.0,-1.0,0.0],7)
                imposed_veclocity = cond.CalculateOnIntegrationPoints(KratosParticle.MPC_IMPOSED_VELOCITY, SomeProcessInfo)
                self.assertVectorAlmostEqual(imposed_veclocity[0],[1.0,-1.0,1.1],7)
                imposed_acceleration = cond.CalculateOnIntegrationPoints(KratosParticle.MPC_IMPOSED_ACCELERATION, SomeProcessInfo)
                self.assertVectorAlmostEqual(imposed_acceleration[0],[1.0,-1.0,2.1],7)
            #base_condition members
            xg = cond.CalculateOnIntegrationPoints(KratosParticle.MPC_COORD, SomeProcessInfo)
            self.assertVectorAlmostEqual(xg[0],[1.5,-1.0,2.1])
            velocity = cond.CalculateOnIntegrationPoints(KratosParticle.MPC_VELOCITY, SomeProcessInfo)
            self.assertVectorAlmostEqual(velocity[0],[1.5,-1.0,2.45])
            acceleration = cond.CalculateOnIntegrationPoints(KratosParticle.MPC_ACCELERATION, SomeProcessInfo)
            self.assertVectorAlmostEqual(acceleration[0],[1.5,-1.12,2.45])

    def _transfer_conditions(self, dimension ):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("mp_dirichlet_conditions")

        rank = data_comm.Rank()
        size = data_comm.Size()
        send_conditions = []

        for i in range(size):
            send_conditions.append( KratosMultiphysics.ConditionsArray() )

        if rank is 0: #Sender
            self._create_condition(mp, dimension, condition_type = "dirichlet", id=1)
            self._create_condition(mp, dimension, condition_type = "dirichlet", id=2)
            for i in range(size):
                if i is not rank:
                    SomeProcessInfo = KratosMultiphysics.ProcessInfo()
                    for cond in mp.Conditions:
                        #Set pseudo-variables
                        self._assign_pseudo_variables(cond, "dirichlet")
                        send_conditions[i].append(cond)
        if rank is 1: #Sender
            self._create_condition(mp, dimension, condition_type = "neumann", id=3)
            for i in range(size):
                if i is not rank:
                    for cond in mp.Conditions:
                        #Set pseudo-variables
                        self._assign_pseudo_variables(cond, "neumann")
                        send_conditions[i].append(cond)

        KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(mp)
        #Send elements from 0 to all other
        KratosParticle.TransferConditions(mp, send_conditions)

        if rank is 0:
            self.assertEqual(mp.NumberOfConditions(),1)
            self._check_conditions(mp, dimension)
        elif rank is 1:
            self.assertEqual(mp.NumberOfConditions(),2)
            self._check_conditions(mp, dimension)
        else:
            self.assertEqual(mp.NumberOfConditions(),3)
            self._check_conditions(mp, dimension)

    def test_transfer_conditions(self):
        self._transfer_conditions(dimension=2)

    #TODO: ADD third dimension

if __name__ == '__main__':
    KratosUnittest.main()