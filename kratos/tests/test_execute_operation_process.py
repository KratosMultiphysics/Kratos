import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.execute_operation_process import ExecuteOperationProcess

class TestExactIntegration(KratosUnittest.TestCase):
    class TestOperation(Kratos.Operation):
        def __init__(self, model: Kratos.Model = None, params: Kratos.Parameters = None):
            super().__init__()

            if model != None and params != None:
                self.model = model
                self.parameters = params

                default_parameters = Kratos.Parameters("""{
                    "model_part_name" : ""
                }""")
                self.parameters.ValidateAndAssignDefaults(default_parameters)
                self.model_part = self.model[self.parameters["model_part_name"].GetString()]
                self.model_part.SetValue(Kratos.PRESSURE, 0)
                self.is_initialized = True
            else:
                self.is_initialized = False

        def Create(self, model: Kratos.Model, params: Kratos.Parameters):
            return TestExactIntegration.TestOperation(model, params)

        def Execute(self):
            if (self.is_initialized):
                self.model_part.SetValue(Kratos.PRESSURE, self.model_part.GetValue(Kratos.PRESSURE) + 1)

    @classmethod
    def setUpClass(cls):
        Kratos.Registry.AddItem("Operations.KratosMultiphysics.TestOperation.Prototype", TestExactIntegration.TestOperation())

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test")

    def test_execute_operation_process(self):
        settings = Kratos.Parameters("""{
            "operation_name"    : "Operations.KratosMultiphysics.TestOperation",
            "operation_settings": {
                "model_part_name" : "test"
            },
            "execution_points"  : [
                "execute_initialize",
                "execute_before_solution_loop",
                "execute_initialize_solution_step",
                "execute_finalize_solution_step",
                "execute_before_output_step",
                "execute_after_output_step",
                "execute_finalize"
            ],
            "echo_level"        : 0
        }""")

        process = ExecuteOperationProcess(self.model, settings)
        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()
        process.ExecuteInitializeSolutionStep()
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteBeforeOutputStep()
        process.ExecuteAfterOutputStep()
        process.ExecuteFinalize()

        self.assertEqual(self.model["test"].GetValue(Kratos.PRESSURE), 7)

    @classmethod
    def tearDownClass(cls):
        Kratos.Registry.RemoveItem("Operations.KratosMultiphysics.TestOperation")


if __name__ == '__main__':
    KratosUnittest.main()
