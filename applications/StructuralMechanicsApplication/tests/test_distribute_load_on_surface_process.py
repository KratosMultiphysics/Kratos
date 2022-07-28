import KratosMultiphysics as KM

import KratosMultiphysics.StructuralMechanicsApplication as KSM
from KratosMultiphysics.StructuralMechanicsApplication.distribute_load_on_surface_process import Factory
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestDistributeLoadOnSurfaceProcess(KratosUnittest.TestCase):

    def test_load_on_surface_distribution(self):
        """Distribute a single force on 2 triangles and 1 quad with a total area of 2.0."""

        current_model = KM.Model()
        mp = current_model.CreateModelPart("main")

        #create nodes
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        mp.CreateNewNode(3, 1.0, 1.0, 0.0)
        mp.CreateNewNode(4, 0.0, 1.0, 0.0)
        mp.CreateNewNode(5, 1.0, 1.0, 1.0)
        mp.CreateNewNode(6, 1.0, 0.0, 1.0)

        #ensure that the property 1 is created
        prop = mp.GetProperties()[1]

        cond1 = mp.CreateNewCondition("SurfaceLoadCondition3D3N", 1, [1,2,3], prop)
        cond2 = mp.CreateNewCondition("SurfaceLoadCondition3D3N", 2, [1,3,4], prop)
        cond3 = mp.CreateNewCondition("SurfaceLoadCondition3D4N", 3, [2,6,5,3], prop)

        settings = KM.Parameters("""{
            "Parameters" : {
                "model_part_name": "main",
                "load": [1.0, 2.0, 3.0]
            }
        }""")

        process = Factory(settings, current_model)

        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()
        process.ExecuteInitializeSolutionStep()

        surface_load_1 = cond1.GetValue(KSM.SURFACE_LOAD)
        surface_load_2 = cond2.GetValue(KSM.SURFACE_LOAD)
        surface_load_3 = cond3.GetValue(KSM.SURFACE_LOAD)
        total_force = KM.Vector([1.0, 2.0, 3.0])

        area1 = cond1.GetGeometry().Area()
        area2 = cond2.GetGeometry().Area()
        area3 = cond3.GetGeometry().Area()
        total_area = area1 + area2 + area3

        self.assertVectorAlmostEqual(surface_load_1, total_force / total_area)
        self.assertVectorAlmostEqual(surface_load_2, total_force / total_area)
        self.assertVectorAlmostEqual(surface_load_3, total_force / total_area)
        self.assertVectorAlmostEqual(surface_load_1 * area1 + surface_load_2 * area2 + surface_load_3 * area3, total_force)

        process.ExecuteFinalizeSolutionStep()
        process.ExecuteBeforeOutputStep()
        process.ExecuteAfterOutputStep()
        process.ExecuteFinalize()


if __name__ == '__main__':
    KratosUnittest.main()
