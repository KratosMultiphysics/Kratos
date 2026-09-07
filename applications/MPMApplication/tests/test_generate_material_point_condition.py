import KratosMultiphysics

import KratosMultiphysics.MPMApplication as KratosMPM
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestGenerateMaterialPointCondition(KratosUnittest.TestCase):

    def _generate_material_point_condition_and_check(self, current_model, dimension, geometry_element, num_material_points, expected_num_material_points):
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
            condition.SetValue(KratosMPM.MATERIAL_POINTS_PER_CONDITION, num_material_points)
            condition.SetValue(KratosMPM.MPC_BOUNDARY_CONDITION_TYPE, 1)
            if (geometry_element == "Point"):
                condition.SetValue(KratosMPM.MPC_IS_NEUMANN, True)
                condition.SetValue(KratosMPM.POINT_LOAD, [0.0,0,0])

        # Create element and nodes for initial meshes
        sub_mp = initial_mesh_model_part.CreateSubModelPart("test")
        sub_mp.GetProperties()[1].SetValue(KratosMPM.MATERIAL_POINTS_PER_ELEMENT, 4)

        # Generate MP Conditions
        KratosMPM.GenerateMaterialPointCondition(grid_model_part, initial_mesh_model_part, material_point_model_part)

        # Check total number of element
        material_point_counter = material_point_model_part.NumberOfConditions()
        self.assertEqual(expected_num_material_points,material_point_counter)

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

    ###########################################################
    ## Point2D - automatic, 1, and default

    def test_GenerateMaterialPointConditionPoint2DAutomatic(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=2, geometry_element="Point", num_material_points=0, expected_num_material_points=1)

    def test_GenerateMaterialPointConditionPoint2D1P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=2, geometry_element="Point", num_material_points=1, expected_num_material_points=1)

    def test_GenerateMaterialPointConditionPoint2DDefault(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=2, geometry_element="Point", num_material_points=50, expected_num_material_points=1)

    ###########################################################
    ## Line2D - 1 and 2, 3, 4, 5, and wrong number

    def test_GenerateMaterialPointConditionLine2D1P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=2, geometry_element="Line", num_material_points=1, expected_num_material_points=1)

    def test_GenerateMaterialPointConditionLine2D2P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=2, geometry_element="Line", num_material_points=2, expected_num_material_points=2)

    def test_GenerateMaterialPointConditionLine2D3P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=2, geometry_element="Line", num_material_points=3, expected_num_material_points=3)

    def test_GenerateMaterialPointConditionLine2D4P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=2, geometry_element="Line", num_material_points=4, expected_num_material_points=4)

    def test_GenerateMaterialPointConditionLine2D5P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=2, geometry_element="Line", num_material_points=5, expected_num_material_points=5)

    # Expected failure: removed default behaviour (which consisted of using 1 material point per line condition)
    @KratosUnittest.expectedFailure
    def test_GenerateMaterialPointConditionLine2DWrongNumber0P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=2, geometry_element="Line", num_material_points=0, expected_num_material_points=1)

    # Expected failure: removed default behaviour (which consisted of using 1 material point per line condition)
    @KratosUnittest.expectedFailure
    def test_GenerateMaterialPointConditionLine2DWrongNumber50P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=2, geometry_element="Line", num_material_points=50, expected_num_material_points=1)

    ###########################################################
    ## Point3D - automatic, 1, and default

    def test_GenerateMaterialPointConditionPoint3DAutomatic(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Point", num_material_points=0, expected_num_material_points=1)

    def test_GenerateMaterialPointConditionPoint3D1P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Point", num_material_points=1, expected_num_material_points=1)

    def test_GenerateMaterialPointConditionPoint3DDefault(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Point", num_material_points=50, expected_num_material_points=1)

    ###########################################################
    ## Line3D - 1, 2, 3, 4, 5, and wrong number

    def test_GenerateMaterialPointConditionLine3D1P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Line", num_material_points=1, expected_num_material_points=1)

    def test_GenerateMaterialPointConditionLine3D2P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Line", num_material_points=2, expected_num_material_points=2)

    def test_GenerateMaterialPointConditionLine3D3P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Line", num_material_points=3, expected_num_material_points=3)

    def test_GenerateMaterialPointConditionLine3D4P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Line", num_material_points=4, expected_num_material_points=4)

    def test_GenerateMaterialPointConditionLine3D5P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Line", num_material_points=5, expected_num_material_points=5)

    # Expected failure: removed default behaviour (which consisted of using 1 material point per line condition)
    @KratosUnittest.expectedFailure
    def test_GenerateMaterialPointConditionLine3DWrongNumber0P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Line", num_material_points=0, expected_num_material_points=1)

    # Expected failure: removed default behaviour (which consisted of using 1 material point per line condition)
    @KratosUnittest.expectedFailure
    def test_GenerateMaterialPointConditionLine3DWrongNumber50P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Line", num_material_points=50, expected_num_material_points=1)

    ###########################################################
    ## Triangle3D - 1, 3, 6, 12, and wrong number

    def test_GenerateMaterialPointConditionTriangle3D1P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Triangle", num_material_points=1, expected_num_material_points=1)

    def test_GenerateMaterialPointConditionTriangle3D3P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Triangle", num_material_points=3, expected_num_material_points=3)

    def test_GenerateMaterialPointConditionTriangle3D6P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Triangle", num_material_points=6, expected_num_material_points=6)

    def test_GenerateMaterialPointConditionTriangle3D12P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Triangle", num_material_points=12, expected_num_material_points=12)

    # Expected failure: removed default behaviour (which consisted of using 1 material point condition per triangle)
    @KratosUnittest.expectedFailure
    def test_GenerateMaterialPointConditionTriangle3DWrongNumber0P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Triangle", num_material_points=0, expected_num_material_points=1)

    # Expected failure: removed default behaviour (which consisted of using 1 material point condition per triangle)
    @KratosUnittest.expectedFailure
    def test_GenerateMaterialPointConditionTriangle3DWrongNumber50P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Triangle", num_material_points=50, expected_num_material_points=1)

    ###########################################################
    ## Quadrilateral3D - 1 ,4, 9, 16 and wrong number

    def test_GenerateMaterialPointConditionQuadrilateral3D4N(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_material_points=4, expected_num_material_points=4)

    def test_GenerateMaterialPointConditionQuadrilateral3D9N(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_material_points=9, expected_num_material_points=9)

    def test_GenerateMaterialPointConditionQuadrilateral3D16N(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_material_points=16, expected_num_material_points=16)

    # Expected failure: removed default behaviour (which consisted of using 1 material point condition per quadrilateral)
    @KratosUnittest.expectedFailure
    def test_GenerateMaterialPointConditionQuadrilateral3DWrongNumber0P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_material_points=0, expected_num_material_points=1)

    # Expected failure: removed default behaviour (which consisted of using 1 material point condition per quadrilateral)
    @KratosUnittest.expectedFailure
    def test_GenerateMaterialPointConditionQuadrilateral3DWrongNumber50P(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition_and_check(current_model, dimension=3, geometry_element="Quadrilateral", num_material_points=50, expected_num_material_points=1)


if __name__ == '__main__':
    KratosUnittest.main()
