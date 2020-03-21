from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestParticleEraseProcess(KratosUnittest.TestCase):

    def _generate_particle_element_and_check(self, current_model):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        dimension = 3

        # Initialize model part
        ## Material model part definition
        material_point_model_part = current_model.CreateModelPart("dummy_name")
        material_point_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)

        ## Initial material model part definition
        initial_mesh_model_part = current_model.CreateModelPart("Initial_dummy_name")
        initial_mesh_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)

        ## Grid model part definition
        grid_model_part = current_model.CreateModelPart("Background_Grid")
        grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)

        # Create element and nodes
        sub_mp = initial_mesh_model_part.CreateSubModelPart("test")
        sub_mp.GetProperties()[1].SetValue(KratosParticle.PARTICLES_PER_ELEMENT, 4)
        self._create_nodes(sub_mp)
        self._create_elements(sub_mp)

        # Create background element and nodes
        background_sub_mp = grid_model_part.CreateSubModelPart("test2")
        self._create_nodes(background_sub_mp)
        self._create_elements(background_sub_mp)
        self._create_conditions(background_sub_mp)

        # Generate MP Elements and Conditions
        KratosParticle.GenerateMaterialPointElement(grid_model_part, initial_mesh_model_part, material_point_model_part, False, False)
        KratosParticle.GenerateMaterialPointCondition(grid_model_part, initial_mesh_model_part, material_point_model_part)

    def _create_nodes(self, initial_mp):
        initial_mp.CreateNewNode(1, -0.5, -0.5, 0.0)
        initial_mp.CreateNewNode(2,  0.5, -0.5, 0.0)
        initial_mp.CreateNewNode(3,  0.5,  0.5, 0.0)
        initial_mp.CreateNewNode(4, -0.5,  0.5, 0.0)
        initial_mp.CreateNewNode(5, -0.5, -0.5, 1.0)
        initial_mp.CreateNewNode(6,  0.5, -0.5, 1.0)
        initial_mp.CreateNewNode(7,  0.5,  0.5, 1.0)
        initial_mp.CreateNewNode(8, -0.5,  0.5, 1.0)

    def _create_elements(self, initial_mp):
        initial_mp.CreateNewElement("UpdatedLagrangian3D8N", 1, [1,2,3,4,5,6,7,8], initial_mp.GetProperties()[1])
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, initial_mp.Elements)

    def _create_conditions(self, initial_mp):
        initial_mp.CreateNewCondition("SurfaceCondition3D4N", 1, [2,4,8,6], initial_mp.GetProperties()[1])
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.BOUNDARY, True, initial_mp.Conditions)
        for condition in initial_mp.Conditions:
            condition.SetValue(KratosParticle.PARTICLES_PER_CONDITION, 0)

    def _search_element(self, current_model):
        # Default
        max_num_results = 1000
        specific_tolerance = 1.e-5

        # Get model part
        material_point_model_part = current_model.GetModelPart("dummy_name")
        grid_model_part           = current_model.GetModelPart("Background_Grid")

        # Search element
        KratosParticle.SearchElement(grid_model_part, material_point_model_part, max_num_results, specific_tolerance)


    def test_ParticleEraseOutsideGivenDomain(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model)

        # Get model part
        material_point_model_part = current_model.GetModelPart("dummy_name")

        # Check initial total number of element
        particle_counter = material_point_model_part.NumberOfElements()
        self.assertEqual(particle_counter, 8)

        # Move particle
        for mpm in material_point_model_part.Elements:
            new_coordinate = mpm.GetValue(KratosParticle.MP_COORD) + [0.3, 0.23, 0.22]
            mpm.SetValue(KratosParticle.MP_COORD, new_coordinate)

        # Check outside given domain
        for mpm in material_point_model_part.Elements:
            new_coordinate = mpm.GetValue(KratosParticle.MP_COORD)
            if(new_coordinate[0] < -0.5 or new_coordinate[0] > 0.5 or new_coordinate[1] < -0.5 or new_coordinate[1] > 0.5 or new_coordinate[2] < -0.5 or new_coordinate[2] > 0.5 ):
                mpm.Set(KratosMultiphysics.TO_ERASE, True)

        # Initiate process
        process = KratosParticle.ParticleEraseProcess(material_point_model_part)

        # Execute
        process.Execute()

        # Check total number of element
        particle_counter = material_point_model_part.NumberOfElements()
        self.assertEqual(particle_counter, 1)
        expected_id = 9
        for mpm in material_point_model_part.Elements:
            self.assertEqual(mpm.Id, expected_id)

    def test_ParticleEraseBySearch(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model)

        # Get model part
        material_point_model_part = current_model.GetModelPart("dummy_name")

        # Check initial total number of element
        particle_counter = material_point_model_part.NumberOfElements()
        self.assertEqual(particle_counter, 8)

        # Move particle
        for mpm in material_point_model_part.Elements:
            new_coordinate = mpm.GetValue(KratosParticle.MP_COORD) + [0.3, 0.23, 0.22]
            mpm.SetValue(KratosParticle.MP_COORD, new_coordinate)

        # Call Search
        self._search_element(current_model)

        # Initiate process
        process = KratosParticle.ParticleEraseProcess(material_point_model_part)

        # Execute
        process.Execute()

        # Check total number of element
        particle_counter = material_point_model_part.NumberOfElements()
        self.assertEqual(particle_counter, 1)
        expected_id = 9
        for mpm in material_point_model_part.Elements:
            self.assertEqual(mpm.Id, expected_id)

    def test_ParticleConditionEraseOutsideGivenDomain(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model)

        # Get model part
        material_point_model_part = current_model.GetModelPart("dummy_name")

        # Check initial number of condition
        particle_counter = material_point_model_part.NumberOfConditions()
        self.assertEqual(particle_counter, 4)

        # Move particle
        for mpc in material_point_model_part.Conditions:
            new_coordinate = mpc.GetValue(KratosParticle.MPC_COORD) + [-0.5, 0.5, 0.5]
            mpc.SetValue(KratosParticle.MPC_COORD, new_coordinate)

        # Check outside given domain
        for mpc in material_point_model_part.Conditions:
            new_coordinate = mpc.GetValue(KratosParticle.MPC_COORD)
            if(new_coordinate[0] < -0.5 or new_coordinate[0] > 0.5 or new_coordinate[1] < -0.5 or new_coordinate[1] > 0.5 or new_coordinate[2] < -0.5 or new_coordinate[2] > 0.5 ):
                mpc.Set(KratosMultiphysics.TO_ERASE, True)

        # Initiate process
        process = KratosParticle.ParticleEraseProcess(material_point_model_part)

        # Execute
        process.Execute()

        # Check total number of condition
        particle_counter = material_point_model_part.NumberOfConditions()
        self.assertEqual(particle_counter, 1)
        expected_id = 11
        for mpc in material_point_model_part.Conditions:
            self.assertEqual(mpc.Id, expected_id)

    def test_ParticleConditionEraseBySearch(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model)

        # Get model part
        material_point_model_part = current_model.GetModelPart("dummy_name")

        # Check initial number of condition
        particle_counter = material_point_model_part.NumberOfConditions()
        self.assertEqual(particle_counter, 4)

        # Move particle
        for mpc in material_point_model_part.Conditions:
            new_coordinate = mpc.GetValue(KratosParticle.MPC_COORD) + [-0.5, 0.5, 0.5]
            mpc.SetValue(KratosParticle.MPC_COORD, new_coordinate)

        # Call Search
        self._search_element(current_model)

        # Initiate process
        process = KratosParticle.ParticleEraseProcess(material_point_model_part)

        # Execute
        process.Execute()

        # Check total number of condition
        particle_counter = material_point_model_part.NumberOfConditions()
        self.assertEqual(particle_counter, 1)
        expected_id = 11
        for mpc in material_point_model_part.Conditions:
            self.assertEqual(mpc.Id, expected_id)


if __name__ == '__main__':
    KratosUnittest.main()