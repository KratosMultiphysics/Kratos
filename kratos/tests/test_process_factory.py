import KratosMultiphysics
from KratosMultiphysics.process_factory import KratosProcessFactory
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestProcessFactory(KratosUnittest.TestCase):
    def test_GetMapOfRequiredHistoricalVariables(self):
        model = KratosMultiphysics.Model()
        processes_list = KratosMultiphysics.Parameters('''[
        {
            "name" : "Processes.KratosMultiphysics.Process",
            "Parameters"    : {}
        },
        {
            "name" : "Processes.KratosMultiphysics.AssignScalarVariableProcess",
            "Parameters"    : {
                "model_part_name" : "main_model_part",
                "variable_name"   : "VELOCITY"
            }
        }
        ]''')
        factory = KratosProcessFactory(model)
        map_variables = factory.GetMapOfRequiredHistoricalVariables(processes_list)
        self.assertEqual(len(map_variables),1)
        for key,variables in map_variables.items():
            self.assertEqual(key, "main_model_part")
            self.assertEqual(len(variables),1)
            for variable in variables:
                self.assertEqual(variable,"VELOCITY")


if __name__ == "__main__":
    KratosUnittest.main()
