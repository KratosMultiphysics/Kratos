from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestGenerateMPMParticle(KratosUnittest.TestCase):

    def _generate_particle_element_and_check(self, current_model, dimension, geometry_element, num_particle, expected_num_particle):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

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
        self._create_nodes(sub_mp, dimension, geometry_element)
        self._create_elements(sub_mp,dimension, geometry_element)

        # Initialize linear_solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()

        # Initialize element
        if geometry_element == "Triangle":
            if (dimension == 2):
                new_element = KratosParticle.CreateUpdatedLagragian2D3N()
            else:
                new_element = KratosParticle.CreateUpdatedLagragian3D4N()
        elif geometry_element == "Quadrilateral":
            if (dimension == 2):
                new_element = KratosParticle.CreateUpdatedLagragian2D4N()
            else:
                new_element = KratosParticle.CreateUpdatedLagragian3D8N()

        # Initialize solver
        if(dimension==2):
            solver = KratosParticle.MPM2D(grid_model_part, initial_material_model_part, material_model_part, linear_solver, new_element, False, "static", geometry_element, num_particle, False, False)
        else:
            solver = KratosParticle.MPM3D(grid_model_part, initial_material_model_part, material_model_part, linear_solver, new_element, False, "static", geometry_element, num_particle, False, False)

        # Check total number of element
        particle_counter = len(material_model_part.Elements)
        self.assertEqual(expected_num_particle,particle_counter)

    def _create_nodes(self, initial_mp, dimension, geometry_element):
        if geometry_element == "Triangle":
            initial_mp.CreateNewNode(1, 0.0, 0.0, 0.0)
            initial_mp.CreateNewNode(2, 1.0, 0.0, 0.0)
            initial_mp.CreateNewNode(3, 0.0, 1.0, 0.0)
            if (dimension == 3):
                initial_mp.CreateNewNode(4, 0.0, 0.0, 1.0)
        elif geometry_element == "Quadrilateral":
            initial_mp.CreateNewNode(1, -0.5, -0.5, 0.0)
            initial_mp.CreateNewNode(2,  0.5, -0.5, 0.0)
            initial_mp.CreateNewNode(3,  0.5,  0.5, 0.0)
            initial_mp.CreateNewNode(4, -0.5,  0.5, 0.0)
            if (dimension == 3):
                initial_mp.CreateNewNode(5, -0.5, -0.5, 1.0)
                initial_mp.CreateNewNode(6,  0.5, -0.5, 1.0)
                initial_mp.CreateNewNode(7,  0.5,  0.5, 1.0)
                initial_mp.CreateNewNode(8, -0.5,  0.5, 1.0)

    def _create_elements(self, initial_mp, dimension, geometry_element):
        if geometry_element == "Triangle":
            if (dimension == 2):
                initial_mp.CreateNewElement("UpdatedLagrangian2D3N", 1, [1,2,3], initial_mp.GetProperties()[1])
            if (dimension == 3):
                initial_mp.CreateNewElement("UpdatedLagrangian3D4N", 1, [1,2,3,4], initial_mp.GetProperties()[1])
        elif geometry_element == "Quadrilateral":
            if (dimension == 2):
                initial_mp.CreateNewElement("UpdatedLagrangian2D4N", 1, [1,2,3,4], initial_mp.GetProperties()[1])
            if (dimension == 3):
                initial_mp.CreateNewElement("UpdatedLagrangian3D8N", 1, [1,2,3,4,5,6,7,8], initial_mp.GetProperties()[1])

        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, initial_mp.Elements)

    def test_GenerateMPMParticleTriangle2D_1(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Triangle", num_particle=1, expected_num_particle=1)

    def test_GenerateMPMParticleTriangle2D_3(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Triangle", num_particle=3, expected_num_particle=3)

    def test_GenerateMPMParticleTriangle2D_6(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Triangle", num_particle=6, expected_num_particle=6)

    def test_GenerateMPMParticleTriangle2D_16(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Triangle", num_particle=16, expected_num_particle=16)

    def test_GenerateMPMParticleTriangle3D_1(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Triangle", num_particle=1, expected_num_particle=1)

    def test_GenerateMPMParticleTriangle3D_3(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Triangle", num_particle=3, expected_num_particle=4)

    def test_GenerateMPMParticleTriangle3D_6(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Triangle", num_particle=6, expected_num_particle=14)

    def test_GenerateMPMParticleTriangle3D_16(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Triangle", num_particle=16, expected_num_particle=16)

    def test_GenerateMPMParticleQuadrilateral2D_1(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Quadrilateral", num_particle=1, expected_num_particle=1)

    def test_GenerateMPMParticleQuadrilateral2D_4(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Quadrilateral", num_particle=4, expected_num_particle=4)

    def test_GenerateMPMParticleQuadrilateral2D_9(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Quadrilateral", num_particle=9, expected_num_particle=9)

    def test_GenerateMPMParticleQuadrilateral2D_16(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Quadrilateral", num_particle=16, expected_num_particle=16)

    def test_GenerateMPMParticleQuadrilateral3D_1(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=1, expected_num_particle=1)

    def test_GenerateMPMParticleQuadrilateral3D_4(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=4, expected_num_particle=8)

    def test_GenerateMPMParticleQuadrilateral3D_9(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=9, expected_num_particle=27)

    def test_GenerateMPMParticleQuadrilateral3D_16(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=16, expected_num_particle=64)

if __name__ == '__main__':
    KratosUnittest.main()
