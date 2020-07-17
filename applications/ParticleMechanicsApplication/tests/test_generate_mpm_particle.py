from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestGenerateMPMParticle(KratosUnittest.TestCase):

    def _generate_particle_element_and_check(self, current_model, dimension, geometry_element, num_particle, expected_num_particle):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

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

        # Create element and nodes for background grids
        sub_background = grid_model_part.CreateSubModelPart("test_background")
        self._create_nodes(sub_background, dimension, geometry_element)
        self._create_elements(sub_background,dimension, geometry_element)

        # Create element and nodes for initial meshes
        sub_mp = initial_mesh_model_part.CreateSubModelPart("test")
        sub_mp.GetProperties()[1].SetValue(KratosParticle.PARTICLES_PER_ELEMENT, num_particle)
        self._create_nodes(sub_mp, dimension, geometry_element)
        self._create_elements(sub_mp,dimension, geometry_element)

        # Generate MP Elements
        KratosParticle.GenerateMaterialPointElement(grid_model_part, initial_mesh_model_part, material_point_model_part, False)

        # Check total number of element
        particle_counter = material_point_model_part.NumberOfElements()
        self.assertEqual(expected_num_particle,particle_counter)

    def _generate_particle_element_and_check_mp_volume(self, current_model, dimension, geometry_element, num_particle, expected_mp_volume):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

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

        # Create element and nodes for background grids
        sub_background = grid_model_part.CreateSubModelPart("test_background")
        self._create_nodes(sub_background, dimension, geometry_element)
        self._create_elements(sub_background,dimension, geometry_element)

        # Create element and nodes for initial meshes
        sub_mp = initial_mesh_model_part.CreateSubModelPart("test")
        sub_mp.GetProperties()[1].SetValue(KratosParticle.PARTICLES_PER_ELEMENT, num_particle)
        self._create_nodes(sub_mp, dimension, geometry_element)
        self._create_elements(sub_mp,dimension, geometry_element)

        # Generate MP Elements
        KratosParticle.GenerateMaterialPointElement(grid_model_part, initial_mesh_model_part, material_point_model_part, False)

        # Check volume of first material point
        for mp in material_point_model_part.Elements:
            mp_volume = mp.CalculateOnIntegrationPoints(KratosParticle.MP_VOLUME, grid_model_part.ProcessInfo)[0]
            self.assertEqual(expected_mp_volume,mp_volume)
            break

    def _create_nodes(self, initial_mp, dimension, geometry_element):
        if geometry_element == "Triangle":
            initial_mp.CreateNewNode(1, 0.0, 0.0, 0.0)
            initial_mp.CreateNewNode(2, 1.0, 0.0, 0.0)
            initial_mp.CreateNewNode(3, 0.0, 1.0, 0.0)
            if (dimension == 3):
                initial_mp.CreateNewNode(4, 0.0, 0.0, 1.0)
        elif geometry_element == "TriangleSkew":
            initial_mp.CreateNewNode(1, 0.0, 0.0, 0.0)
            initial_mp.CreateNewNode(2, 2.0, 0.0, 0.0)
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
        elif geometry_element == "QuadrilateralSkew":
            initial_mp.CreateNewNode(1, -0.5, -0.5, 0.0)
            initial_mp.CreateNewNode(2,  1.5, -0.5, 0.0)
            initial_mp.CreateNewNode(3,  0.5,  0.5, 0.0)
            initial_mp.CreateNewNode(4, -0.5,  0.5, 0.0)
            if (dimension == 3):
                initial_mp.CreateNewNode(5, -0.5, -0.5, 1.0)
                initial_mp.CreateNewNode(6,  0.5, -0.5, 1.0)
                initial_mp.CreateNewNode(7,  0.5,  0.5, 1.0)
                initial_mp.CreateNewNode(8, -0.5,  0.5, 1.0)

    def _create_elements(self, initial_mp, dimension, geometry_element):
        if geometry_element == "Triangle" or geometry_element == "TriangleSkew":
            if (dimension == 2):
                initial_mp.CreateNewElement("Element2D3N", 1, [1,2,3], initial_mp.GetProperties()[1])
            if (dimension == 3):
                initial_mp.CreateNewElement("Element3D4N", 1, [1,2,3,4], initial_mp.GetProperties()[1])
        elif geometry_element == "Quadrilateral" or geometry_element == "QuadrilateralSkew":
            if (dimension == 2):
                initial_mp.CreateNewElement("Element2D4N", 1, [1,2,3,4], initial_mp.GetProperties()[1])
            if (dimension == 3):
                initial_mp.CreateNewElement("Element3D8N", 1, [1,2,3,4,5,6,7,8], initial_mp.GetProperties()[1])

        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, initial_mp.Elements)

    def test_GenerateMPMParticleTriangle2D1P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Triangle", num_particle=1, expected_num_particle=1)

    def test_GenerateMPMParticleTriangle2D3P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Triangle", num_particle=3, expected_num_particle=3)

    def test_GenerateMPMParticleTriangle2D6P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Triangle", num_particle=6, expected_num_particle=6)

    def test_GenerateMPMParticleTriangle2D12P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Triangle", num_particle=12, expected_num_particle=12)

    def test_GenerateMPMParticleTriangle2D16P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Triangle", num_particle=16, expected_num_particle=16)

    def test_GenerateMPMParticleTriangle2D33P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Triangle", num_particle=33, expected_num_particle=33)

    def test_GenerateMPMParticleTriangle2DDefault(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Triangle", num_particle=50, expected_num_particle=3)

    def test_GenerateMPMParticleTriangle3D1P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Triangle", num_particle=1, expected_num_particle=1)

    def test_GenerateMPMParticleTriangle3D4P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Triangle", num_particle=3, expected_num_particle=4)

    def test_GenerateMPMParticleTriangle3D14P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Triangle", num_particle=6, expected_num_particle=14)

    def test_GenerateMPMParticleTriangle3D24P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Triangle", num_particle=12, expected_num_particle=24)

    def test_GenerateMPMParticleTriangle3DDefault(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Triangle", num_particle=50, expected_num_particle=4)

    def test_GenerateMPMParticleQuadrilateral2D1P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Quadrilateral", num_particle=1, expected_num_particle=1)

    def test_GenerateMPMParticleQuadrilateral2D4P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Quadrilateral", num_particle=4, expected_num_particle=4)

    def test_GenerateMPMParticleQuadrilateral2D9P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Quadrilateral", num_particle=9, expected_num_particle=9)

    def test_GenerateMPMParticleQuadrilateral2D16P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Quadrilateral", num_particle=16, expected_num_particle=16)

    def test_GenerateMPMParticleQuadrilateral2DDefault(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=2, geometry_element="Quadrilateral", num_particle=50, expected_num_particle=4)

    def test_GenerateMPMParticleQuadrilateral3D1P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=1, expected_num_particle=1)

    def test_GenerateMPMParticleQuadrilateral3D8P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=4, expected_num_particle=8)

    def test_GenerateMPMParticleQuadrilateral3D27P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=9, expected_num_particle=27)

    def test_GenerateMPMParticleQuadrilateral3D64P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=16, expected_num_particle=64)

    def test_GenerateMPMParticleQuadrilateral3DDefault(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=50, expected_num_particle=8)

    # Tests for the correct computation of material point volume in the material point generator
    def test_GenerateMPMParticleQuadrilateral2DSkew(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check_mp_volume(current_model, dimension=2, geometry_element="QuadrilateralSkew", num_particle=4, expected_mp_volume=0.44716878364870316)

    def test_GenerateMPMParticleQuadrilateral3DSkew(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check_mp_volume(current_model, dimension=3, geometry_element="QuadrilateralSkew", num_particle=4, expected_mp_volume=0.20275105849101815)
    
    def test_GenerateMPMParticleTriangle2DSkew(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check_mp_volume(current_model, dimension=2, geometry_element="TriangleSkew", num_particle=3, expected_mp_volume=0.3333333333333333)
    
    def test_GenerateMPMParticleTriangle3DSkew(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check_mp_volume(current_model, dimension=3, geometry_element="TriangleSkew", num_particle=3, expected_mp_volume=0.08333333333333333)

    def test_GenerateMPMParticleQuadrilateral2DNotSkew(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check_mp_volume(current_model, dimension=2, geometry_element="Quadrilateral", num_particle=4, expected_mp_volume=0.25)

    def test_GenerateMPMParticleQuadrilateral3DNotSkew(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check_mp_volume(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=4, expected_mp_volume=0.12499999999999993)
    
    def test_GenerateMPMParticleTriangle2DNotSkew(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check_mp_volume(current_model, dimension=2, geometry_element="Triangle", num_particle=3, expected_mp_volume=0.16666666666666666)
    
    def test_GenerateMPMParticleTriangle3DNotSkew(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_element_and_check_mp_volume(current_model, dimension=3, geometry_element="Triangle", num_particle=3, expected_mp_volume=0.041666666666666664)

if __name__ == '__main__':
    KratosUnittest.main()
