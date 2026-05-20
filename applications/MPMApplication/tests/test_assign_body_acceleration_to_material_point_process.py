import math

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.MPMApplication as KratosMPM
from KratosMultiphysics.MPMApplication.assign_body_acceleration_to_material_point_process import (
    AssignBodyAccelerationToMaterialPointProcess
)


class TestAssignBodyAccelerationToMaterialPointProcess(KratosUnittest.TestCase):

    def test_assigns_time_and_space_dependent_value_to_volume_acceleration(self):
        model = KratosMultiphysics.Model()
        mp_model_part = self._CreateMaterialPointModelPart(model)
        mp_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.5)

        element = next(iter(mp_model_part.Elements))
        element.SetValuesOnIntegrationPoints(
            KratosMPM.MP_ACCELERATION,
            [[9.0, 8.0, 7.0]],
            mp_model_part.ProcessInfo)

        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "MPM_Material",
            "value"           : ["2.0*t", "sin(t)", "x+y+t"]
        }""")

        process = AssignBodyAccelerationToMaterialPointProcess(model, settings)
        process.ExecuteInitializeSolutionStep()

        mp_coord = element.CalculateOnIntegrationPoints(KratosMPM.MP_COORD, mp_model_part.ProcessInfo)[0]
        volume_acceleration = element.CalculateOnIntegrationPoints(
            KratosMPM.MP_VOLUME_ACCELERATION,
            mp_model_part.ProcessInfo)[0]
        acceleration = element.CalculateOnIntegrationPoints(
            KratosMPM.MP_ACCELERATION,
            mp_model_part.ProcessInfo)[0]

        self.assertVectorAlmostEqual(
            volume_acceleration,
            [1.0, math.sin(0.5), mp_coord[0] + mp_coord[1] + 0.5])
        self.assertVectorAlmostEqual(acceleration, [9.0, 8.0, 7.0])

    def test_keeps_legacy_modulus_component_settings(self):
        model = KratosMultiphysics.Model()
        mp_model_part = self._CreateMaterialPointModelPart(model)
        mp_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.25)

        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "MPM_Material",
            "modulus"         : "cos(t)",
            "component"       : [1.0, "2.0*t", 0.0],
            "set_initial_mp_acceleration": true
        }""")

        process = AssignBodyAccelerationToMaterialPointProcess(model, settings)
        process.ExecuteBeforeSolutionLoop()

        element = next(iter(mp_model_part.Elements))
        expected_value = [math.cos(0.25), 0.5 * math.cos(0.25), 0.0]
        volume_acceleration = element.CalculateOnIntegrationPoints(
            KratosMPM.MP_VOLUME_ACCELERATION,
            mp_model_part.ProcessInfo)[0]
        acceleration = element.CalculateOnIntegrationPoints(
            KratosMPM.MP_ACCELERATION,
            mp_model_part.ProcessInfo)[0]

        self.assertVectorAlmostEqual(volume_acceleration, expected_value)
        self.assertVectorAlmostEqual(acceleration, expected_value)

    def _CreateMaterialPointModelPart(self, model):
        grid_model_part = model.CreateModelPart("Background_Grid")
        grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        initial_model_part = model.CreateModelPart("Initial_MPM_Material")
        initial_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        mp_model_part = model.CreateModelPart("MPM_Material")
        mp_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        self._CreateQuadrilateral(grid_model_part.CreateSubModelPart("Grid"))
        material_sub_model_part = initial_model_part.CreateSubModelPart("Material")
        self._CreateQuadrilateral(material_sub_model_part)

        properties = material_sub_model_part.GetProperties()[1]
        properties.SetValue(KratosMPM.MATERIAL_POINTS_PER_ELEMENT, 1)
        properties.SetValue(KratosMultiphysics.DENSITY, 1000.0)
        properties.SetValue(KratosMultiphysics.THICKNESS, 1.0)

        KratosMPM.GenerateMaterialPointElement(grid_model_part, initial_model_part, mp_model_part, False)
        return mp_model_part

    @staticmethod
    def _CreateQuadrilateral(model_part):
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        model_part.CreateNewNode(4, 0.0, 1.0, 0.0)
        model_part.CreateNewElement("Element2D4N", 1, [1, 2, 3, 4], model_part.GetProperties()[1])
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, model_part.Elements)


if __name__ == '__main__':
    KratosUnittest.main()
