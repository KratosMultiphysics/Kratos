import KratosMultiphysics
from KratosMultiphysics import KratosUnittest
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.MPMApplication as KratosMPM
data_comm = KratosMultiphysics.DataCommunicator.GetDefault()

class TestTransferConditions(KratosUnittest.TestCase):
    ''' This class provides all required methods to test the MPM_MPI_Utilities::TransferConditions function.
        Tests in 2D and 3D with dirichlet, neumann and dirichlet coupling interface conditons are performed.
        New developed conditions must be added here.
    '''
    def _create_nodes(self, mp, dimension):
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,0.0,1.0,0.0)
        mp.CreateNewNode(3,1.0,0.0,0.0)
        if (dimension == 3):
            mp.CreateNewNode(4, 0.0, 0.0, 1.0)

    def _create_particle_condition(self, mp, dimension, condition_type, condition_id):
        #Create nodes
        self._create_nodes(mp, dimension)

        if dimension == 2:
            if condition_type == "dirichlet":
                mp.CreateNewCondition("MPMParticlePenaltyDirichletCondition2D3N", condition_id, [1, 2, 3], mp.GetProperties()[1])
            if condition_type == "neumann":
                mp.CreateNewCondition("MPMParticlePointLoadCondition2D3N", condition_id, [1, 2, 3], mp.GetProperties()[1])
        if dimension == 3:
            if condition_type == "dirichlet":
                mp.CreateNewCondition("MPMParticlePenaltyDirichletCondition3D4N", condition_id, [1, 2, 3, 4], mp.GetProperties()[1])
            if condition_type == "neumann":
                mp.CreateNewCondition("MPMParticlePointLoadCondition3D4N", condition_id, [1, 2, 3, 4], mp.GetProperties()[1])

    def _assign_pseudo_variables(self, cond, condition_type):
        process_info = KratosMultiphysics.ProcessInfo()
        normal_vector = [KratosMultiphysics.Vector([1.0,0.0,0.0])]
        cond.SetValuesOnIntegrationPoints(KratosMPM.MPC_NORMAL,normal_vector,process_info)
        xg = [KratosMultiphysics.Vector([1.5,-1.0,2.1])]
        cond.SetValuesOnIntegrationPoints(KratosMPM.MPC_COORD,xg,process_info)
        velocity = [KratosMultiphysics.Vector([1.5,-1.0,2.45])]
        cond.SetValuesOnIntegrationPoints(KratosMPM.MPC_VELOCITY,velocity,process_info)
        acceleration = [KratosMultiphysics.Vector([1.5,-1.12,2.45])]
        cond.SetValuesOnIntegrationPoints(KratosMPM.MPC_ACCELERATION,acceleration,process_info)
        if condition_type == "dirichlet" or condition_type == "coupling":
            cond.SetValuesOnIntegrationPoints(KratosMPM.PENALTY_FACTOR,[100.0],process_info)
            displacement = [KratosMultiphysics.Vector([1.22,-1.11,0.0])]
            cond.SetValuesOnIntegrationPoints(KratosMPM.MPC_DISPLACEMENT,displacement,process_info)
            imposed_displacement = [KratosMultiphysics.Vector([1.0,-1.0,0.0])]
            cond.SetValuesOnIntegrationPoints(KratosMPM.MPC_IMPOSED_DISPLACEMENT,imposed_displacement,process_info)
            imposed_velocity = [KratosMultiphysics.Vector([1.0,-1.0,1.1])]
            cond.SetValuesOnIntegrationPoints(KratosMPM.MPC_IMPOSED_VELOCITY,imposed_velocity,process_info)
            imposed_acceleration = [KratosMultiphysics.Vector([1.0,-1.0,2.1])]
            cond.SetValuesOnIntegrationPoints(KratosMPM.MPC_IMPOSED_ACCELERATION,imposed_acceleration,process_info)
        else:
            point_load = [KratosMultiphysics.Vector([3.3,4.4,5.5])]
            cond.SetValuesOnIntegrationPoints(KratosMPM.POINT_LOAD,point_load ,process_info)

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
            if dimension == 2:
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
            unit_normal = cond.CalculateOnIntegrationPoints(KratosMPM.MPC_NORMAL, process_info)
            self.assertVectorAlmostEqual(unit_normal[0],[1.0,0.0,0.0],7)
            if(cond.Info() == "Condition #3"):
                #point_load_condition members
                point_load = cond.CalculateOnIntegrationPoints(KratosMPM.POINT_LOAD, process_info)
                self.assertVectorAlmostEqual(point_load[0],[3.3,4.4,5.5])
            else:
                #penalty_dirichlet_condition members
                penalty_factor = cond.CalculateOnIntegrationPoints(KratosMPM.PENALTY_FACTOR, process_info)
                self.assertAlmostEqual(penalty_factor[0], 100.0, 7)
                #base_dirichlet_condition members
                displacement = cond.CalculateOnIntegrationPoints(KratosMPM.MPC_DISPLACEMENT, process_info)
                self.assertVectorAlmostEqual(displacement[0],[1.22,-1.11,0.0],7)
                imposed_displacement = cond.CalculateOnIntegrationPoints(KratosMPM.MPC_IMPOSED_DISPLACEMENT, process_info)
                self.assertVectorAlmostEqual(imposed_displacement[0],[1.0,-1.0,0.0],7)
                imposed_veclocity = cond.CalculateOnIntegrationPoints(KratosMPM.MPC_IMPOSED_VELOCITY, process_info)
                self.assertVectorAlmostEqual(imposed_veclocity[0],[1.0,-1.0,1.1],7)
                imposed_acceleration = cond.CalculateOnIntegrationPoints(KratosMPM.MPC_IMPOSED_ACCELERATION, process_info)
                self.assertVectorAlmostEqual(imposed_acceleration[0],[1.0,-1.0,2.1],7)
            #base_condition members
            xg = cond.CalculateOnIntegrationPoints(KratosMPM.MPC_COORD, process_info)
            self.assertVectorAlmostEqual(xg[0],[1.5,-1.0,2.1])
            velocity = cond.CalculateOnIntegrationPoints(KratosMPM.MPC_VELOCITY, process_info)
            self.assertVectorAlmostEqual(velocity[0],[1.5,-1.0,2.45])
            acceleration = cond.CalculateOnIntegrationPoints(KratosMPM.MPC_ACCELERATION, process_info)
            self.assertVectorAlmostEqual(acceleration[0],[1.5,-1.12,2.45])

    def _transfer_conditions(self, dimension, condition_type_2 ):
        ''' Two dirichlet conditions are created in rank=0 and send to all other processes.
            One neumann or coupled/interface condition is created in rank=1 and send to all other processes.
            The test is passed if all receiving processes hold the correct conditions. '''
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("mp_dirichlet_conditions")

        rank = data_comm.Rank()
        size = data_comm.Size()
        send_conditions = []

        for i in range(size):
            send_conditions.append( KratosMultiphysics.ConditionsArray() )

        if rank == 0: #Sender
            self._create_particle_condition(mp, dimension, condition_type="dirichlet", condition_id=1)
            self._create_particle_condition(mp, dimension, condition_type="dirichlet", condition_id=2)
            for i in range(size):
                if i != rank:
                    for cond in mp.Conditions:
                        #Set pseudo-variables
                        self._assign_pseudo_variables(cond, "dirichlet")
                        send_conditions[i].append(cond)
        if rank == 1: #Sender
            if condition_type_2 == "neumann":
                condition_id = 3
            else:
                condition_id = 4
            self._create_particle_condition(mp, dimension, condition_type_2, condition_id)
            for i in range(size):
                if i != rank:
                    for cond in mp.Conditions:
                        #Set pseudo-variables
                        self._assign_pseudo_variables(cond, condition_type_2)
                        send_conditions[i].append(cond)

        KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(mp)
        # Exchange elements
        KratosMPM.MPM_MPI_Utilities.TransferConditions(mp, send_conditions)

        # Check
        if rank == 0:
            self.assertEqual(mp.NumberOfConditions(),1)
            self._check_conditions(mp, dimension)
        elif rank == 1:
            self.assertEqual(mp.NumberOfConditions(),2)
            self._check_conditions(mp, dimension)
        else:
            self.assertEqual(mp.NumberOfConditions(),3)
            self._check_conditions(mp, dimension)

    def test_transfer_conditions2D_dirichlet_neumann(self):
        self._transfer_conditions(dimension=2, condition_type_2="neumann")

    def test_transfer_conditions3D_dirichlet_neumann(self):
        self._transfer_conditions(dimension=3, condition_type_2="neumann")

if __name__ == '__main__':
    KratosUnittest.main()