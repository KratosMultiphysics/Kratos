import KratosMultiphysics

import KratosMultiphysics.MPMApplication as KratosMPM
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestSearchMaterialPointCondition(KratosUnittest.TestCase):

    def _generate_material_point_condition(self, current_model, dimension, geometry_element):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # Initialize model part
        ## Material model part definition
        material_point_model_part = current_model.CreateModelPart("dummy_name")
        material_point_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)
        self.process_info = material_point_model_part.ProcessInfo

        ## Initial material model part definition
        initial_mesh_model_part = current_model.CreateModelPart("Initial_dummy_name")
        initial_mesh_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)

        ## Grid model part definition
        grid_model_part = current_model.CreateModelPart("Background_Grid")
        grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, dimension)

        # Create Background Grid
        sub_background = grid_model_part.CreateSubModelPart("test")
        self._create_background_nodes(sub_background, dimension, geometry_element)

        self._create_background_elements(sub_background,dimension, geometry_element)

        self._create_nodes(sub_background, dimension)

        self._create_conditions(sub_background,dimension)
        for condition in grid_model_part.Conditions:
            condition.SetValue(KratosMPM.MATERIAL_POINTS_PER_CONDITION, 1)
            condition.SetValue(KratosMPM.MPC_IS_NEUMANN, True)
            condition.SetValue(KratosMPM.POINT_LOAD, [1.0, 0.0, 0.0])

        # Set active
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, initial_mesh_model_part.Elements)

        # Generate MP Conditions
        KratosMPM.GenerateMaterialPointCondition(grid_model_part, initial_mesh_model_part, material_point_model_part)


    def _create_nodes(self, model_part, dimension):
        if (dimension == 2):
            model_part.CreateNewNode(13, 0.1, 0.1, 0.0)
        if (dimension == 3):
            model_part.CreateNewNode(13, 0.1, 0.1, 0.1)


    def _create_background_nodes(self, model_part, dimension, geometry_element):
        if geometry_element == "Triangle":
            model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
            model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
            model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
            model_part.CreateNewNode(5, 1.0, 1.0, 0.0)
            if (dimension == 3):
                model_part.CreateNewNode(4, 0.0, 0.0, 1.0)
                model_part.CreateNewNode(6, 1.0, 0.0, 1.0)
                model_part.CreateNewNode(7, 0.0, 1.0, 1.0)
                model_part.CreateNewNode(8, 1.0, 1.0, 1.0)
        elif geometry_element == "Quadrilateral":
            model_part.CreateNewNode(1, -0.5, -0.5, 0.0)
            model_part.CreateNewNode(2,  0.5, -0.5, 0.0)
            model_part.CreateNewNode(3,  0.5,  0.5, 0.0)
            model_part.CreateNewNode(4, -0.5,  0.5, 0.0)
            model_part.CreateNewNode(9 , 1.5, -0.5, 0.0)
            model_part.CreateNewNode(10, 1.5,  0.5, 0.0)

            if (dimension == 3):
                model_part.CreateNewNode(5, -0.5, -0.5, 1.0)
                model_part.CreateNewNode(6,  0.5, -0.5, 1.0)
                model_part.CreateNewNode(7,  0.5,  0.5, 1.0)
                model_part.CreateNewNode(8, -0.5,  0.5, 1.0)
                model_part.CreateNewNode(11, 1.5, -0.5, 1.0)
                model_part.CreateNewNode(12, 1.5,  0.5, 1.0)


    def _create_conditions(self, model_part, dimension):
        if (dimension == 2):
            model_part.CreateNewCondition("PointCondition2D1N", 1, [13], model_part.GetProperties()[1])
        else:
            model_part.CreateNewCondition("PointCondition3D1N", 1, [13], model_part.GetProperties()[1])

        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.BOUNDARY, True, model_part.Conditions)


    def _create_background_elements(self, model_part, dimension, geometry_element):
        if geometry_element == "Triangle":
            if (dimension == 2):
                model_part.CreateNewElement("Element2D3N", 1, [1,2,3], model_part.GetProperties()[1])
                model_part.CreateNewElement("Element2D3N", 2, [2,3,5], model_part.GetProperties()[1])
            if (dimension == 3):
                model_part.CreateNewElement("Element3D4N", 1, [1,2,3,4], model_part.GetProperties()[1])
                model_part.CreateNewElement("Element3D4N", 2, [2,8,4,6], model_part.GetProperties()[1])
                model_part.CreateNewElement("Element3D4N", 3, [4,8,3,7], model_part.GetProperties()[1])
                model_part.CreateNewElement("Element3D4N", 4, [2,5,3,8], model_part.GetProperties()[1])
                model_part.CreateNewElement("Element3D4N", 5, [8,3,2,4], model_part.GetProperties()[1])
        elif geometry_element == "Quadrilateral":
            if (dimension == 2):
                model_part.CreateNewElement("Element2D4N", 1, [1,2,3,4], model_part.GetProperties()[1])
                model_part.CreateNewElement("Element2D4N", 2, [2,9,10,3], model_part.GetProperties()[1])
            if (dimension == 3):
                model_part.CreateNewElement("Element3D8N", 1, [1,2,3,4,5,6,7,8], model_part.GetProperties()[1])
                model_part.CreateNewElement("Element3D8N", 2, [2,9,10,3,6,11,12,7], model_part.GetProperties()[1])


    def _move_and_search_condition(self, current_model, new_coordinate, max_num_results = 1000, specific_tolerance = 1.e-5):
        # Get model part
        material_point_model_part = current_model.GetModelPart("dummy_name")
        grid_model_part           = current_model.GetModelPart("Background_Grid")

        # Apply before search
        for mpc in material_point_model_part.Conditions:
            mpc.SetValuesOnIntegrationPoints(KratosMPM.MPC_COORD, [new_coordinate], self.process_info)

        # Search element
        KratosMPM.SearchElement(grid_model_part, material_point_model_part, max_num_results, specific_tolerance)

    def _check_connectivity(self, current_model, expected_connectivity_node=[]):
        # Get model part
        material_point_model_part = current_model.GetModelPart("dummy_name")
        grid_model_part           = current_model.GetModelPart("Background_Grid")

        # Check the searched node as expected connectivity
        if not expected_connectivity_node:
            for mpc in material_point_model_part.Conditions:
                self.assertEqual(mpc.GetNodes(), [])
        else:
            for mpc in material_point_model_part.Conditions:
                if (mpc.GetGeometry().PointsNumber() == 0):
                    print("Search Element: Condition was not found!")
                    for i in range (len(expected_connectivity_node)):
                        self.assertEqual(mpc.GetNode(i).Id, ())
                else:
                    for i in range (len(expected_connectivity_node)):
                        self.assertEqual(mpc.GetNode(i).Id, grid_model_part.GetNode(expected_connectivity_node[i]).Id)
                        self.assertEqual(mpc.GetNode(i).X, grid_model_part.GetNode(expected_connectivity_node[i]).X)
                        self.assertEqual(mpc.GetNode(i).Y, grid_model_part.GetNode(expected_connectivity_node[i]).Y)
                        self.assertEqual(mpc.GetNode(i).Z, grid_model_part.GetNode(expected_connectivity_node[i]).Z)

    def test_SearchMaterialPointConditionTriangle2D(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition(current_model, dimension=2, geometry_element="Triangle")

        new_coordinate = [0.25, 0.25, 0.0]
        self._move_and_search_condition(current_model, new_coordinate)
        self._check_connectivity(current_model, [1,2,3])

        new_coordinate = [0.50001, 0.50001, 0.0]
        self._move_and_search_condition(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,3,5])

        new_coordinate = [1.00001, 1.00001, 0.0]
        self._move_and_search_condition(current_model, new_coordinate)
        self._check_connectivity(current_model)

    def test_SearchMaterialPointConditionTriangle3D(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition(current_model, dimension=3, geometry_element="Triangle")

        new_coordinate = [0.5, 0.25, 0.20]
        self._move_and_search_condition(current_model, new_coordinate)
        self._check_connectivity(current_model, [1,2,3,4])

        new_coordinate = [0.90, 0.55, 0.90]
        self._move_and_search_condition(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,8,4,6])

        new_coordinate = [0.10, 0.90, 0.55]
        self._move_and_search_condition(current_model, new_coordinate)
        self._check_connectivity(current_model, [4,8,3,7])

        new_coordinate = [0.90, 0.90, 0.55]
        self._move_and_search_condition(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,5,3,8])

        new_coordinate = [0.50, 0.50, 0.50]
        self._move_and_search_condition(current_model, new_coordinate)
        self._check_connectivity(current_model, [8,3,2,4])

        new_coordinate = [1.0001, 1.0001, 1.0001]
        self._move_and_search_condition(current_model, new_coordinate)
        self._check_connectivity(current_model)

    def test_SearchMaterialPointConditionQuadrilateral2D(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition(current_model, dimension=2, geometry_element="Quadrilateral")

        new_coordinate = [-0.11111, 0.12345, 0.0]
        self._move_and_search_condition(current_model, new_coordinate)
        self._check_connectivity(current_model, [1,2,3,4])

        new_coordinate = [0.6, 0.12345, 0.0]
        self._move_and_search_condition(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,9,10,3])

        new_coordinate = [1.00001, 1.00001, 0.0]
        self._move_and_search_condition(current_model, new_coordinate)
        self._check_connectivity(current_model)

    def test_SearchMaterialPointConditionQuadrilateral3D(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_condition(current_model, dimension=3, geometry_element="Quadrilateral")

        new_coordinate = [0.5, 0.25, 0.20]
        self._move_and_search_condition(current_model, new_coordinate)
        self._check_connectivity(current_model, [1,2,3,4,5,6,7,8])

        new_coordinate = [0.7, 0.35, 0.3]
        self._move_and_search_condition(current_model, new_coordinate)
        self._check_connectivity(current_model, [2,9,10,3,6,11,12,7])

        new_coordinate = [0.50001, 0.50001, 0.50001]
        self._move_and_search_condition(current_model, new_coordinate)
        self._check_connectivity(current_model)


if __name__ == '__main__':
    KratosUnittest.main()
