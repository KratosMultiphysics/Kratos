import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

class TestPredictor(KratosUnittest.TestCase):
    def test_predictor(self):
        parameter_file_name = "test_parameters.json"
        cs_data_structure = ImportDataStructure(parameter_file_name)

        m = 10
        dz = 3.0
        a0 = 1.0
        a1 = 2.0
        a2 = 3.0
        a3 = 4.0
        a4 = 5.0
        interface_settings = cs_data_structure.Parameters('{"wall": "AREA"}')

        # Create interface
        variable = vars(KM)["AREA"]
        model = cs_data_structure.Model()
        model_part = model.CreateModelPart("wall")
        model_part.AddNodalSolutionStepVariable(variable)
        for i in range(m):
            model_part.CreateNewNode(i, 0.0, 0.0, i * dz)
        step = 0
        for node in model_part.Nodes:
            node.SetSolutionStepValue(variable, step, a0)
        interface = CoSimulationInterface(model, interface_settings)

        # Create predictor
        parameter_file_name = "predictors/test_cubic.json"
        with open(parameter_file_name, 'r') as parameter_file:
            settings = cs_data_structure.Parameters(parameter_file.read())
        predictor_cubic = cs_tools.CreateInstance(settings)
        predictor_cubic.Initialize(interface)
        interface_as_array = interface.GetNumpyArray()

        # Test predictor: a linear relation should be predicted in the same way
        # by linear, quadratic and cubic predictors
        predictor_cubic.InitializeSolutionStep()
        interface.SetNumpyArray(a1 * interface_as_array)
        predictor_cubic.Update(interface)
        predictor_cubic.FinalizeSolutionStep()
        predictor_cubic.InitializeSolutionStep()
        interface.SetNumpyArray(a2 * interface_as_array)
        predictor_cubic.Update(interface)
        predictor_cubic.FinalizeSolutionStep()
        predictor_cubic.InitializeSolutionStep()
        interface.SetNumpyArray(a3 * interface_as_array)
        predictor_cubic.Update(interface)
        predictor_cubic.FinalizeSolutionStep()

        predictor_cubic.InitializeSolutionStep()
        prediction_linear = predictor_cubic.Linear(interface).GetNumpyArray()
        prediction_quadratic = predictor_cubic.Quadratic(interface).GetNumpyArray()
        prediction_cubic = predictor_cubic.Cubic(interface).GetNumpyArray()
        for i in range(m):
            self.assertAlmostEqual(a4, prediction_linear[i])
            self.assertAlmostEqual(a4, prediction_quadratic[i])
            self.assertAlmostEqual(a4, prediction_cubic[i])

        # Test predictor: error if no update
        with self.assertRaises(Exception):
            predictor_cubic.InitializeSolutionStep()
            predictor_cubic.FinalizeSolutionStep()

        # Test predictor: error if updated twice
        with self.assertRaises(Exception):
            predictor_cubic.InitializeSolutionStep()
            prediction = predictor_cubic.Predict(interface)
            prediction = predictor_cubic.Predict(interface)
            predictor_cubic.FinalizeSolutionStep()

        # Test predictor: error if prediction after update
        with self.assertRaises(Exception):
            predictor_cubic.InitializeSolutionStep()
            prediction = predictor_cubic.Update(interface)
            prediction = predictor_cubic.Predict(interface)
            predictor_cubic.FinalizeSolutionStep()


if __name__ == '__main__':
    KratosUnittest.main()
