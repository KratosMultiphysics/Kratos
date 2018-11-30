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
        material_model_part = current_model.CreateModelPart("dummy_name")
        material_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)

        ## Initial material model part definition
        initial_material_model_part = current_model.CreateModelPart("Initial_dummy_name")
        initial_material_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)

        ## Grid model part definition
        grid_model_part = current_model.CreateModelPart("Background_Grid")
        grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)

        # Create element and nodes
        sub_mp = initial_material_model_part.CreateSubModelPart("test")
        self._create_nodes(sub_mp)
        self._create_elements(sub_mp)

        # Create background element and nodes
        background_sub_mp = grid_model_part.CreateSubModelPart("test")
        self._create_nodes(background_sub_mp)
        self._create_elements(background_sub_mp)

        # Initialize linear_solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()

        # Initialize element
        new_element = KratosParticle.CreateUpdatedLagragian3D8N()

        # Initialize solver
        self.solver = KratosParticle.MPM3D(grid_model_part, initial_material_model_part, material_model_part, linear_solver, new_element, False, "static", "Quadrilateral", 4, False, False)

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

    def _search_element(self, current_model):
        # Default
        max_num_results = 1000
        specific_tolerance = 1.e-5

        # Get model part
        material_model_part = current_model.GetModelPart("dummy_name")
        grid_model_part     = current_model.GetModelPart("Background_Grid")

        # Search element
        self.solver.SearchElement(grid_model_part, material_model_part, max_num_results, specific_tolerance)


    def test_ParticleEraseOutsideGivenDomain(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model)

        # Get model part
        material_model_part = current_model.GetModelPart("dummy_name")

        # Move particle
        for mpm in material_model_part.Elements:
            new_coordinate = mpm.GetValue(KratosParticle.MP_COORD) + [0.3, 0.23, 0.22]
            mpm.SetValue(KratosParticle.MP_COORD, new_coordinate)

        # Check outside given domain
        for mpm in material_model_part.Elements:
            new_coordinate = mpm.GetValue(KratosParticle.MP_COORD)
            if(new_coordinate[0] < -0.5 or new_coordinate[0] > 0.5 or new_coordinate[1] < -0.5 or new_coordinate[1] > 0.5 or new_coordinate[2] < -0.5 or new_coordinate[2] > 0.5 ):
                mpm.Set(KratosMultiphysics.TO_ERASE, True)

        # Initiate process
        process = KratosParticle.ParticleEraseProcess(material_model_part)

        # Execute
        process.Execute()

        # Check total number of element
        particle_counter = len(material_model_part.Elements)
        self.assertEqual(particle_counter, 1)
        expected_id = 9
        for mpm in material_model_part.Elements:
            self.assertEqual(mpm.Id, expected_id)

    def test_ParticleEraseBySearch(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model)

        # Get model part
        material_model_part = current_model.GetModelPart("dummy_name")

        # Move particle
        for mpm in material_model_part.Elements:
            new_coordinate = mpm.GetValue(KratosParticle.MP_COORD) + [0.3, 0.23, 0.22]
            mpm.SetValue(KratosParticle.MP_COORD, new_coordinate)

        # Call Search
        self._search_element(current_model)

        # Initiate process
        process = KratosParticle.ParticleEraseProcess(material_model_part)

        # Execute
        process.Execute()

        # Check total number of element
        particle_counter = len(material_model_part.Elements)
        self.assertEqual(particle_counter, 1)
        expected_id = 9
        for mpm in material_model_part.Elements:
            self.assertEqual(mpm.Id, expected_id)

if __name__ == '__main__':
    KratosUnittest.main()