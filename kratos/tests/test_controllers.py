import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestControllers(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = KM.Model()
        cls.model_part = cls.model.CreateModelPart("test")

    def setUp(self) -> None:
        self.model_part.ProcessInfo[KM.STEP] = 0
        self.model_part.ProcessInfo[KM.TIME] = 0

    def test_TemporalControllerStep(self):
        params = KM.Parameters("""{
            "model_part_name"     : "test",
            "output_control_type" : "step",
            "output_interval"     : 3
        }""")

        # testing from start
        controller = KM.OutputController(self.model, params)
        for _ in range(100):
            self.model_part.ProcessInfo[KM.STEP] += 1
            self.assertEqual(controller.Evaluate(), self.model_part.ProcessInfo[KM.STEP] % 3 == 0)
            controller.Update()

        # testing for restart
        self.model_part.ProcessInfo[KM.STEP] = 99
        controller = KM.OutputController(self.model, params)
        for _ in range(100):
            self.model_part.ProcessInfo[KM.STEP] += 1
            self.assertEqual(controller.Evaluate(), self.model_part.ProcessInfo[KM.STEP] % 3 == 0)
            controller.Update()

    def test_TemporalControllerTime(self):
        params = KM.Parameters("""{
            "model_part_name"     : "test",
            "output_control_type" : "time",
            "output_interval"     : 3.0
        }""")

        # testing from start
        controller = KM.OutputController(self.model, params)
        for _ in range(100):
            self.model_part.ProcessInfo[KM.TIME] += 1
            self.assertEqual(controller.Evaluate(), self.model_part.ProcessInfo[KM.TIME] % 3 == 0)
            controller.Update()

        # testing for restart
        self.model_part.ProcessInfo[KM.TIME] = 99
        controller = KM.OutputController(self.model, params)
        for _ in range(100):
            self.model_part.ProcessInfo[KM.TIME] += 1
            self.assertEqual(controller.Evaluate(), self.model_part.ProcessInfo[KM.TIME] % 3 == 0)
            controller.Update()

if __name__ == "__main__":
    KratosUnittest.main()
