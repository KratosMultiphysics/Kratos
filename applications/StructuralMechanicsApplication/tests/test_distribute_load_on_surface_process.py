import KratosMultiphysics as KM

import KratosMultiphysics.StructuralMechanicsApplication as KSM
from KratosMultiphysics.StructuralMechanicsApplication.distribute_load_on_surface_process import Factory
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestDistributeLoadOnSurfaceProcess(KratosUnittest.TestCase):

    def test_load_on_surface_distribution(self):
        """distribute on 2 triangles and 1 quad with a total area of 2.0."""

        current_model = KM.Model()
        mp = current_model.CreateModelPart("main")

        # Create nodes
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        mp.CreateNewNode(3, 1.0, 1.0, 0.0)
        mp.CreateNewNode(4, 0.0, 1.0, 0.0)
        mp.CreateNewNode(5, 1.0, 1.0, 1.0)
        mp.CreateNewNode(6, 1.0, 0.0, 1.0)

        # Ensure that the property 1 is created
        prop = mp.GetProperties()[1]

        # Create conditions
        cond1 = mp.CreateNewCondition("SurfaceLoadCondition3D3N", 1, [1,2,3], prop)
        cond2 = mp.CreateNewCondition("SurfaceLoadCondition3D3N", 2, [1,3,4], prop)
        cond3 = mp.CreateNewCondition("SurfaceLoadCondition3D4N", 3, [2,6,5,3], prop)

        # Generate process
        settings = KM.Parameters("""{
            "Parameters" : {
                "model_part_name" : "main",
                "load"            : [1.0, 2.0, 3.0]
            }
        }""")

        process = Factory(settings, current_model)

        # Initial calls of process
        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()
        process.ExecuteInitializeSolutionStep()

        # Getting loads
        surface_load_1 = cond1.GetValue(KSM.SURFACE_LOAD)
        surface_load_2 = cond2.GetValue(KSM.SURFACE_LOAD)
        surface_load_3 = cond3.GetValue(KSM.SURFACE_LOAD)
        total_load = KM.Vector([1.0, 2.0, 3.0])

        # Get total area
        area1 = cond1.GetGeometry().Area()
        area2 = cond2.GetGeometry().Area()
        area3 = cond3.GetGeometry().Area()
        total_area = area1 + area2 + area3
        self.assertAlmostEqual(total_area, 2.0)

        # Compute distributed load
        distributed_load = total_load / total_area
        self.assertVectorAlmostEqual(surface_load_1, distributed_load)
        self.assertVectorAlmostEqual(surface_load_2, distributed_load)
        self.assertVectorAlmostEqual(surface_load_3, distributed_load)
        self.assertVectorAlmostEqual(surface_load_1 * area1 + surface_load_2 * area2 + surface_load_3 * area3, total_load)

        # More calls (does nothing, but checks that it does not break the process)
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteBeforeOutputStep()
        process.ExecuteAfterOutputStep()
        process.ExecuteFinalize()

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()