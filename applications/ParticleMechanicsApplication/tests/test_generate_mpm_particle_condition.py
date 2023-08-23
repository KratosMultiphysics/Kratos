import KratosMultiphysics

import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestGenerateMPMParticleCondition(KratosUnittest.TestCase):

    def _generate_particle_condition_and_check(self, current_model, dimension, geometry_element, num_particle, expected_num_particle):
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
        self._create_condition(sub_background,dimension, geometry_element)
        for condition in grid_model_part.Conditions:
            condition.SetValue(KratosParticle.PARTICLES_PER_CONDITION, num_particle)
            condition.SetValue(KratosParticle.MPC_BOUNDARY_CONDITION_TYPE, 1)


        # Create element and nodes for initial meshes
        sub_mp = initial_mesh_model_part.CreateSubModelPart("test")
        sub_mp.GetProperties()[1].SetValue(KratosParticle.PARTICLES_PER_ELEMENT, 4)

        # Generate MP Conditions
        KratosParticle.GenerateMaterialPointCondition(grid_model_part, initial_mesh_model_part, material_point_model_part)

        # Check total number of element
        particle_counter = material_point_model_part.NumberOfConditions()
        self.assertEqual(expected_num_particle,particle_counter)

    def _create_nodes(self, initial_mp, dimension, geometry_element):
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
        if (dimension == 2):
            initial_mp.CreateNewElement("MPMUpdatedLagrangian2D4N", 1, [1,2,3,4], initial_mp.GetProperties()[1])
        else:
            initial_mp.CreateNewElement("MPMUpdatedLagrangian3D8N", 1, [1,2,3,4,5,6,7,8], initial_mp.GetProperties()[1])

        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, initial_mp.Elements)

    def _create_condition(self, initial_mp, dimension, geometry_element):
        if (dimension == 2):
            if (geometry_element == "Point"):
                initial_mp.CreateNewCondition("PointCondition2D1N", 1, [1], initial_mp.GetProperties()[1])
            elif (geometry_element == "Line"):
                initial_mp.CreateNewCondition("LineCondition2D2N", 1, [1,2], initial_mp.GetProperties()[1])
        else:
            if (geometry_element == "Point"):
                initial_mp.CreateNewCondition("PointCondition3D1N", 1, [1], initial_mp.GetProperties()[1])
            elif (geometry_element == "Line"):
                initial_mp.CreateNewCondition("LineCondition3D2N", 1, [1,2], initial_mp.GetProperties()[1])
            elif (geometry_element == "Triangle"):
                initial_mp.CreateNewCondition("SurfaceCondition3D3N", 1, [1,6,8], initial_mp.GetProperties()[1])
            elif (geometry_element == "Quadrilateral"):
                initial_mp.CreateNewCondition("SurfaceCondition3D4N", 1, [2,4,8,6], initial_mp.GetProperties()[1])

        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.BOUNDARY, True, initial_mp.Conditions)

    ## Point2D - automatic, 1, and default
    def test_GenerateMPMParticleConditionPoint2DAutomatic(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=2, geometry_element="Point", num_particle=0, expected_num_particle=1)

    def test_GenerateMPMParticleConditionPoint2D1P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=2, geometry_element="Point", num_particle=1, expected_num_particle=1)

    def test_GenerateMPMParticleConditionPoint2DDefault(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=2, geometry_element="Point", num_particle=50, expected_num_particle=1)

    ## Line2D - automatic and 2, 3, 4, 5, and default
    def test_GenerateMPMParticleConditionLine2DAutomatic(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=2, geometry_element="Line", num_particle=0, expected_num_particle=1)

    def test_GenerateMPMParticleConditionLine2D1P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=2, geometry_element="Line", num_particle=1, expected_num_particle=1)

    def test_GenerateMPMParticleConditionLine2D2P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=2, geometry_element="Line", num_particle=2, expected_num_particle=2)

    def test_GenerateMPMParticleConditionLine2D3P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=2, geometry_element="Line", num_particle=3, expected_num_particle=3)

    def test_GenerateMPMParticleConditionLine2D4P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=2, geometry_element="Line", num_particle=4, expected_num_particle=4)

    def test_GenerateMPMParticleConditionLine2D5P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=2, geometry_element="Line", num_particle=5, expected_num_particle=5)

    def test_GenerateMPMParticleConditionLine2DDefault(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=2, geometry_element="Line", num_particle=50, expected_num_particle=1)

    ## Point3D - automatic, 1, and default
    def test_GenerateMPMParticleConditionPoint3DAutomatic(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Point", num_particle=0, expected_num_particle=1)

    def test_GenerateMPMParticleConditionPoint3D1P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Point", num_particle=1, expected_num_particle=1)

    def test_GenerateMPMParticleConditionPoint3DDefault(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Point", num_particle=50, expected_num_particle=1)

    ## Line3D - automatic and 2, 3, 4, 5, and default
    def test_GenerateMPMParticleConditionLine3DAutomatic(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Line", num_particle=0, expected_num_particle=1)

    def test_GenerateMPMParticleConditionLine3D1P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Line", num_particle=1, expected_num_particle=1)

    def test_GenerateMPMParticleConditionLine3D2P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Line", num_particle=2, expected_num_particle=2)

    def test_GenerateMPMParticleConditionLine3D3P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Line", num_particle=3, expected_num_particle=3)

    def test_GenerateMPMParticleConditionLine3D4P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Line", num_particle=4, expected_num_particle=4)

    def test_GenerateMPMParticleConditionLine3D5P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Line", num_particle=5, expected_num_particle=5)

    def test_GenerateMPMParticleConditionLine3DDefault(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Line", num_particle=50, expected_num_particle=1)

    ## Triangle3D - automatic, 1, 3, 6, 12, and default
    def test_GenerateMPMParticleConditionTriangle3DAutomatic(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Triangle", num_particle=0, expected_num_particle=1)

    def test_GenerateMPMParticleConditionTriangle3D1P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Triangle", num_particle=1, expected_num_particle=1)

    def test_GenerateMPMParticleConditionTriangle3D3P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Triangle", num_particle=3, expected_num_particle=3)

    def test_GenerateMPMParticleConditionTriangle3D6P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Triangle", num_particle=6, expected_num_particle=6)

    def test_GenerateMPMParticleConditionTriangle3D12P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Triangle", num_particle=12, expected_num_particle=12)

    def test_GenerateMPMParticleConditionTriangle3DDefault(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Triangle", num_particle=50, expected_num_particle=1)

    ## Quadrilateral3D - automatic, 1 ,4, 9, 16 and default
    def test_GenerateMPMParticleConditionQuadrilateral3DAutomatic(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=0, expected_num_particle=1)

    def test_GenerateMPMParticleConditionQuadrilateral3D4N(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=4, expected_num_particle=4)

    def test_GenerateMPMParticleConditionQuadrilateral3D9N(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=9, expected_num_particle=9)

    def test_GenerateMPMParticleConditionQuadrilateral3D16N(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=16, expected_num_particle=16)

    def test_GenerateMPMParticleConditionQuadrilateral3DDefault(self):
        current_model = KratosMultiphysics.Model()
        self._generate_particle_condition_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_particle=50, expected_num_particle=1)


if __name__ == '__main__':
    KratosUnittest.main()
