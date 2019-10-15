import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

class TestPredictorLinear(KratosUnittest.TestCase):
    def test_predictor_linear(self):
        parameter_file_name = "test_parameters.json"
        cs_data_structure = ImportDataStructure(parameter_file_name)

        m = 10
        dz = 2.0
        a0 = 1.0
        a1 = 1.5
        a2 = 2.0
        a3 = 2.5
        interface_settings = cs_data_structure.Parameters('{"wall": "AREA"}')

        # Create interface
        variable = vars(KM)["AREA"]
        model = cs_data_structure.Model()
        model_part = model.CreateModelPart("wall")
        model_part.AddNodalSolutionStepVariable(variable)
        for i in range(m):
            model_part.CreateNewNode(i, 0.0, 0.0, i*dz)
        step = 0
        for node in model_part.Nodes:
            node.SetSolutionStepValue(variable, step, a0)
        interface = CoSimulationInterface(model, interface_settings)

        # Create predictor
        parameter_file_name = "predictors/test_linear.json"
        with open(parameter_file_name, 'r') as parameter_file:
            settings = cs_data_structure.Parameters(parameter_file.read())

        predictor_linear = cs_tools.CreateInstance(settings)
        predictor_linear.Initialize(interface)

        # Test predictor: first prediction needs to be equal to initialized value
        predictor_linear.InitializeSolutionStep()
        prediction = predictor_linear.Predict(interface)
        self.assertIsInstance(prediction, CoSimulationInterface)
        prediction_as_array = prediction.GetNumpyArray()
        for i in range(m):
            self.assertAlmostEqual(a0, prediction_as_array[i])
        interface_as_array = a1*prediction_as_array
        interface.SetNumpyArray(interface_as_array)
        predictor_linear.Update(interface)
        predictor_linear.FinalizeSolutionStep()

        # Test predictor: second prediction needs to be linear
        predictor_linear.InitializeSolutionStep()
        prediction = predictor_linear.Predict(interface)
        self.assertIsInstance(prediction, CoSimulationInterface)
        prediction_as_array = prediction.GetNumpyArray()
        for i in range(m):
            self.assertAlmostEqual(a2, prediction_as_array[i])
        interface_as_array = a3 * prediction_as_array
        interface.SetNumpyArray(interface_as_array)
        predictor_linear.Update(interface)
        predictor_linear.FinalizeSolutionStep()

        # Test predictor: error if no update
        with self.assertRaises(Exception):
            predictor_linear.InitializeSolutionStep()
            predictor_linear.FinalizeSolutionStep()

        # Test predictor: error if updated twice
        with self.assertRaises(Exception):
            predictor_linear.InitializeSolutionStep()
            prediction = predictor_linear.Predict(interface)
            prediction = predictor_linear.Predict(interface)
            predictor_linear.FinalizeSolutionStep()

        # Test predictor: error if prediction after update
        with self.assertRaises(Exception):
            predictor_linear.InitializeSolutionStep()
            prediction = predictor_linear.Update(interface)
            prediction = predictor_linear.Predict(interface)
            predictor_linear.FinalizeSolutionStep()


if __name__ == '__main__':
    KratosUnittest.main()
