import KratosMultiphysics
from KratosMultiphysics import KratosUnittest
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
data_comm = KratosMultiphysics.DataCommunicator.GetDefault()

class TestTransferConditions(KratosUnittest.TestCase):
    ''' This class provides all required methods to test the MPM_MPI_Utilities::TransferConditions function.
        Tests in 2D and 3D with dirichlet and neumann particle conditons are performed.
    '''
    def _create_nodes(self, mp, dimension):
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,0.0,1.0,0.0)
        mp.CreateNewNode(3,1.0,0.0,0.0)
        if (dimension == 3):
            mp.CreateNewNode(4, 0.0, 0.0, 1.0)

    def _create_particle_condition(self, mp, dimension, condition_type, id):
        #Create nodes
        self._create_nodes(mp, dimension)

        # Ensure that the property 1 is created
        mp.GetProperties()[1]

        if dimension == 2:
            if condition_type == "dirichlet":
                cond = mp.CreateNewCondition("MPMParticlePenaltyDirichletCondition2D3N", id, [1, 2, 3], mp.GetProperties()[1])
            if condition_type == "neumann":
                cond = mp.CreateNewCondition("MPMParticlePointLoadCondition2D3N", id, [1, 2, 3], mp.GetProperties()[1])
        if dimension == 3:
            if condition_type == "dirichlet":
                cond = mp.CreateNewCondition("MPMParticlePenaltyDirichletCondition3D4N", id, [1, 2, 3, 4], mp.GetProperties()[1])
            if condition_type == "neumann":
                cond = mp.CreateNewCondition("MPMParticlePointLoadCondition3D4N", id, [1, 2, 3, 4], mp.GetProperties()[1])

    def _assign_pseudo_variables(self, cond, condition_type):
        process_info = KratosMultiphysics.ProcessInfo()
        normal_vector = [KratosMultiphysics.Vector([1.0,0.0,0.0])]
        cond.SetValuesOnIntegrationPoints(KratosParticle.MPC_NORMAL,normal_vector,process_info)
        xg = [KratosMultiphysics.Vector([1.5,-1.0,2.1])]
        cond.SetValuesOnIntegrationPoints(KratosParticle.MPC_COORD,xg,process_info)
        velocity = [KratosMultiphysics.Vector([1.5,-1.0,2.45])]
        cond.SetValuesOnIntegrationPoints(KratosParticle.MPC_VELOCITY,velocity,process_info)
        acceleration = [KratosMultiphysics.Vector([1.5,-1.12,2.45])]
        cond.SetValuesOnIntegrationPoints(KratosParticle.MPC_ACCELERATION,acceleration,process_info)
        if condition_type == "dirichlet":
            cond.SetValuesOnIntegrationPoints(KratosParticle.PENALTY_FACTOR,[100.0],process_info)
            displacement = [KratosMultiphysics.Vector([1.22,-1.11,0.0])]
            cond.SetValuesOnIntegrationPoints(KratosParticle.MPC_DISPLACEMENT,displacement,process_info)
            imposed_displacement = [KratosMultiphysics.Vector([1.0,-1.0,0.0])]
            cond.SetValuesOnIntegrationPoints(KratosParticle.MPC_IMPOSED_DISPLACEMENT,imposed_displacement,process_info)
            imposed_velocity = [KratosMultiphysics.Vector([1.0,-1.0,1.1])]
            cond.SetValuesOnIntegrationPoints(KratosParticle.MPC_IMPOSED_VELOCITY,imposed_velocity,process_info)
            imposed_acceleration = [KratosMultiphysics.Vector([1.0,-1.0,2.1])]
            cond.SetValuesOnIntegrationPoints(KratosParticle.MPC_IMPOSED_ACCELERATION,imposed_acceleration,process_info)
        else:
            point_load = [KratosMultiphysics.Vector([3.3,4.4,5.5])]
            cond.SetValuesOnIntegrationPoints(KratosParticle.POINT_LOAD,point_load ,process_info)

    def _check_conditions(self, mp, dimension):
        process_info = KratosMultiphysics.ProcessInfo()
        for cond in mp.Conditions:
            #Check geometry
            self.assertEqual( cond.GetGeometry().WorkingSpaceDimension(), dimension)
            self.assertEqual( cond.GetGeometry().LocalSpaceDimension(), dimension)
            self.assertEqual( cond.GetGeometry().IntegrationPointsNumber(), 1)
            jacobian = cond.GetGeometry().Jacobian(0)
            shape_functions_values = cond.GetGeometry().ShapeFunctionsValues()
            shape_functions_derivatives = cond.GetGeometry().ShapeFunctionDerivatives(1,0)
            center = cond.GetGeometry().Center()
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
                center_ref = [1.0/3.0, 1.0/3.0, 0.0]
            else:
                jacobian_ref = KratosMultiphysics.Matrix(3,3,0.0)
                jacobian_ref[0,1] = 1.0
                jacobian_ref[1,0] = 1.0
                jacobian_ref[2,2] = 1.0
                shape_functions_values_ref = KratosMultiphysics.Matrix(1,4,0.25)
                shape_functions_derivatives_ref = KratosMultiphysics.Matrix(4,3,0.0)
                shape_functions_derivatives_ref[0,0] = -1.0
                shape_functions_derivatives_ref[0,1] = -1.0
                shape_functions_derivatives_ref[0,2] = -1.0
                shape_functions_derivatives_ref[1,0] = 1.0
                shape_functions_derivatives_ref[2,1] = 1.0
                shape_functions_derivatives_ref[3,2] = 1.0
                center_ref = [0.25,0.25,0.25]
            self.assertMatrixAlmostEqual(jacobian, jacobian_ref, 7)
            self.assertMatrixAlmostEqual(shape_functions_values, shape_functions_values_ref, 7)
            self.assertMatrixAlmostEqual(shape_functions_derivatives, shape_functions_derivatives_ref, 7)
            self.assertVectorAlmostEqual(center, center_ref)
            ##Check condition properties
            unit_normal = cond.CalculateOnIntegrationPoints(KratosParticle.MPC_NORMAL, process_info)
            self.assertVectorAlmostEqual(unit_normal[0],[1.0,0.0,0.0],7)
            if(cond.Info() == "Condition #3"):
                #point_load_condition members
                point_load = cond.CalculateOnIntegrationPoints(KratosParticle.POINT_LOAD, process_info)
                self.assertVectorAlmostEqual(point_load[0],[3.3,4.4,5.5])
            else:
                #penalty_dirichlet_condition members
                penalty_factor = cond.CalculateOnIntegrationPoints(KratosParticle.PENALTY_FACTOR, process_info)
                self.assertAlmostEqual(penalty_factor[0], 100.0, 7)
                #base_dirichlet_condition members
                displacement = cond.CalculateOnIntegrationPoints(KratosParticle.MPC_DISPLACEMENT, process_info)
                self.assertVectorAlmostEqual(displacement[0],[1.22,-1.11,0.0],7)
                imposed_displacement = cond.CalculateOnIntegrationPoints(KratosParticle.MPC_IMPOSED_DISPLACEMENT, process_info)
                self.assertVectorAlmostEqual(imposed_displacement[0],[1.0,-1.0,0.0],7)
                imposed_veclocity = cond.CalculateOnIntegrationPoints(KratosParticle.MPC_IMPOSED_VELOCITY, process_info)
                self.assertVectorAlmostEqual(imposed_veclocity[0],[1.0,-1.0,1.1],7)
                imposed_acceleration = cond.CalculateOnIntegrationPoints(KratosParticle.MPC_IMPOSED_ACCELERATION, process_info)
                self.assertVectorAlmostEqual(imposed_acceleration[0],[1.0,-1.0,2.1],7)
            #base_condition members
            xg = cond.CalculateOnIntegrationPoints(KratosParticle.MPC_COORD, process_info)
            self.assertVectorAlmostEqual(xg[0],[1.5,-1.0,2.1])
            velocity = cond.CalculateOnIntegrationPoints(KratosParticle.MPC_VELOCITY, process_info)
            self.assertVectorAlmostEqual(velocity[0],[1.5,-1.0,2.45])
            acceleration = cond.CalculateOnIntegrationPoints(KratosParticle.MPC_ACCELERATION, process_info)
            self.assertVectorAlmostEqual(acceleration[0],[1.5,-1.12,2.45])

    def _transfer_conditions(self, dimension ):
        ''' Two dirichlet conditions are created in rank=0 and send to all other processes.
            One neumann condition is created in rank=1 and send to all other processesself.
            The test is passed if all receiving processes hold the correct conditions '''
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("mp_dirichlet_conditions")

        rank = data_comm.Rank()
        size = data_comm.Size()
        send_conditions = []

        for i in range(size):
            send_conditions.append( KratosMultiphysics.ConditionsArray() )

        if rank is 0: #Sender
            self._create_particle_condition(mp, dimension, condition_type = "dirichlet", id=1)
            self._create_particle_condition(mp, dimension, condition_type = "dirichlet", id=2)
            for i in range(size):
                if i is not rank:
                    process_info = KratosMultiphysics.ProcessInfo()
                    for cond in mp.Conditions:
                        #Set pseudo-variables
                        self._assign_pseudo_variables(cond, "dirichlet")
                        send_conditions[i].append(cond)
        if rank is 1: #Sender
            self._create_particle_condition(mp, dimension, condition_type = "neumann", id=3)
            for i in range(size):
                if i is not rank:
                    for cond in mp.Conditions:
                        #Set pseudo-variables
                        self._assign_pseudo_variables(cond, "neumann")
                        send_conditions[i].append(cond)

        KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(mp)
        # Exchange elements
        KratosParticle.TransferConditions(mp, send_conditions)

        # Check
        if rank is 0:
            self.assertEqual(mp.NumberOfConditions(),1)
            self._check_conditions(mp, dimension)
        elif rank is 1:
            self.assertEqual(mp.NumberOfConditions(),2)
            self._check_conditions(mp, dimension)
        else:
            self.assertEqual(mp.NumberOfConditions(),3)
            self._check_conditions(mp, dimension)

    def test_transfer_conditions2D(self):
        self._transfer_conditions(dimension=2)

    def test_transfer_conditions3D(self):
        self._transfer_conditions(dimension=3)

if __name__ == '__main__':
    KratosUnittest.main()