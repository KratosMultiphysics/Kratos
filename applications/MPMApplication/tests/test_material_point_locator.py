import KratosMultiphysics

import KratosMultiphysics.MPMApplication as KratosMPM
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestMaterialPointLocator(KratosUnittest.TestCase):

    def _generate_material_point_elements_and_conditions(self, current_model, dimension):
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

        # Create Background Grid nodes, elements and conditions
        grid = grid_model_part.CreateSubModelPart("grid")
        self._create_nodes(grid, dimension)
        self._create_elements(grid, dimension)
        boundary = grid_model_part.CreateSubModelPart("boundary")
        self._create_conditions(boundary, dimension)

        # Create initial mesh element and nodes
        sub_mp = initial_mesh_model_part.CreateSubModelPart("test")
        sub_mp.GetProperties()[1].SetValue(KratosMPM.MATERIAL_POINTS_PER_ELEMENT, 4)
        self._create_nodes(sub_mp, dimension)
        self._create_elements(sub_mp,dimension)

        # Set active
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, initial_mesh_model_part.Elements)
        # Generate MP Elements
        KratosMPM.GenerateMaterialPointElement(grid_model_part, initial_mesh_model_part, material_point_model_part, False)
        # Generate MP Conditions
        KratosMPM.GenerateMaterialPointCondition(grid_model_part, initial_mesh_model_part, material_point_model_part)


    def _create_nodes(self, model_part, dimension):
        model_part.CreateNewNode(1, -0.5, -0.5, 0.0)
        model_part.CreateNewNode(2,  0.5, -0.5, 0.0)
        model_part.CreateNewNode(3,  0.5,  0.5, 0.0)
        model_part.CreateNewNode(4, -0.5,  0.5, 0.0)
        if (dimension == 3):
            model_part.CreateNewNode(5, -0.5, -0.5, 1.0)
            model_part.CreateNewNode(6,  0.5, -0.5, 1.0)
            model_part.CreateNewNode(7,  0.5,  0.5, 1.0)
            model_part.CreateNewNode(8, -0.5,  0.5, 1.0)


    def _create_elements(self, model_part, dimension):
        if (dimension == 2):
            model_part.CreateNewElement("Element2D4N", 1, [1,2,3,4], model_part.GetProperties()[1])
        if (dimension == 3):
            model_part.CreateNewElement("Element3D8N", 1, [1,2,3,4,5,6,7,8], model_part.GetProperties()[1])


    def _create_conditions(self, model_part, dimension):
        if dimension == 2:
            model_part.CreateNewNode( 9, 0.4, -0.5, 0.0)
            model_part.CreateNewNode(10, 0.4,  0.0, 0.0)
            model_part.CreateNewNode(11, 0.4, +0.5, 0.0)
            model_part.CreateNewCondition("LineCondition2D2N", 1, [ 9, 10], model_part.GetProperties()[1])
            model_part.CreateNewCondition("LineCondition2D2N", 2, [10, 11], model_part.GetProperties()[1])
        elif dimension == 3:
            model_part.CreateNewNode( 9, 0.4, -0.5, 0.0)
            model_part.CreateNewNode(10, 0.4,  0.0, 0.0)
            model_part.CreateNewNode(11, 0.4, +0.5, 0.0)
            model_part.CreateNewNode(12, 0.4, -0.5, 1.0)
            model_part.CreateNewNode(13, 0.4,  0.0, 1.0)
            model_part.CreateNewNode(14, 0.4, +0.5, 1.0)
            model_part.CreateNewCondition("SurfaceCondition3D4N", 1, [ 9, 10, 13, 12], model_part.GetProperties()[1])
            model_part.CreateNewCondition("SurfaceCondition3D4N", 2, [10, 11, 14, 13], model_part.GetProperties()[1])
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.BOUNDARY, True, model_part.Conditions)
        for condition in model_part.Conditions:
            condition.SetValue(KratosMPM.MATERIAL_POINTS_PER_CONDITION, 1)
            condition.SetValue(KratosMPM.MPC_BOUNDARY_CONDITION_TYPE, 1)

    def test_FindMaterialPointElement2D(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_elements_and_conditions(current_model, dimension=2)

        model_part = current_model.GetModelPart("dummy_name")

        point = KratosMultiphysics.Point(0.0, 0.0, 0.0)
        entity_id = KratosMPM.BruteForceMaterialPointLocator(model_part).FindElement(point, 1e-1)
        self.assertEqual(entity_id, -1)

        point = KratosMultiphysics.Point(-0.288675, -0.288675, 0.0)
        entity_id = KratosMPM.BruteForceMaterialPointLocator(model_part).FindElement(point, 1e-4)
        self.assertEqual(entity_id, 8)

        point = KratosMultiphysics.Point(-0.288675, -0.288675, 0.0)
        entity_id = KratosMPM.BruteForceMaterialPointLocator(model_part).FindElement(point, 1)
        self.assertEqual(entity_id, 8)

        point = KratosMultiphysics.Point(-0.01, -0.01, 0.0)
        entity_id = KratosMPM.BruteForceMaterialPointLocator(model_part).FindElement(point, 0.4)
        self.assertEqual(entity_id, 8)


    def test_FindMaterialPointElement3D(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_elements_and_conditions(current_model, dimension=3)

        model_part = current_model.GetModelPart("dummy_name")

        point = KratosMultiphysics.Point(0.0, 0.0, 0.0)
        entity_id = KratosMPM.BruteForceMaterialPointLocator(model_part).FindElement(point, 1e-4)
        self.assertEqual(entity_id, -1)

        point = KratosMultiphysics.Point(-0.288675, 0.288675, 0.788675)
        entity_id = KratosMPM.BruteForceMaterialPointLocator(model_part).FindElement(point, 1e-4)
        self.assertEqual(entity_id, 22)

        point = KratosMultiphysics.Point(-0.288675, 0.288675, 0.788675)
        entity_id = KratosMPM.BruteForceMaterialPointLocator(model_part).FindElement(point, 1)
        self.assertEqual(entity_id, 22)

        point = KratosMultiphysics.Point(-0.01, 0.01, 0.51)
        entity_id = KratosMPM.BruteForceMaterialPointLocator(model_part).FindElement(point, 0.5)
        self.assertEqual(entity_id, 22)


    def test_FindMaterialPointCondition2D(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_elements_and_conditions(current_model, dimension=2)

        model_part = current_model.GetModelPart("dummy_name")

        point = KratosMultiphysics.Point(0.0, 0.0, 0.0)
        entity_id = KratosMPM.BruteForceMaterialPointLocator(model_part).FindCondition(point, 1e-1)
        self.assertEqual(entity_id, -1)

        point = KratosMultiphysics.Point(0.4, 0.25, 0.0)
        entity_id = KratosMPM.BruteForceMaterialPointLocator(model_part).FindCondition(point, 1e-4)
        self.assertEqual(entity_id, 9)

        point = KratosMultiphysics.Point(0.4, 0.25, 0.0)
        entity_id = KratosMPM.BruteForceMaterialPointLocator(model_part).FindCondition(point, 1)
        self.assertEqual(entity_id, 9)

        point = KratosMultiphysics.Point(0.4, 0.01, 0.0)
        entity_id = KratosMPM.BruteForceMaterialPointLocator(model_part).FindCondition(point, 0.25)
        self.assertEqual(entity_id, 9)


    def test_FindMaterialPointCondition3D(self):
        current_model = KratosMultiphysics.Model()
        self._generate_material_point_elements_and_conditions(current_model, dimension=3)

        model_part = current_model.GetModelPart("dummy_name")

        point = KratosMultiphysics.Point(0.0, 0.0, 0.0)
        entity_id = KratosMPM.BruteForceMaterialPointLocator(model_part).FindCondition(point, 1e-4)
        self.assertEqual(entity_id, -1)

        point = KratosMultiphysics.Point(0.4, 0.25, 0.5)
        entity_id = KratosMPM.BruteForceMaterialPointLocator(model_part).FindCondition(point, 1e-4)
        self.assertEqual(entity_id, 16)

        point = KratosMultiphysics.Point(0.4, 0.25, 0.5)
        entity_id = KratosMPM.BruteForceMaterialPointLocator(model_part).FindCondition(point, 1)
        self.assertEqual(entity_id, 16)

        point = KratosMultiphysics.Point(0.4, 0.25, 0.5)
        entity_id = KratosMPM.BruteForceMaterialPointLocator(model_part).FindCondition(point, 0.5)
        self.assertEqual(entity_id, 16)


if __name__ == '__main__':
    KratosUnittest.main()
