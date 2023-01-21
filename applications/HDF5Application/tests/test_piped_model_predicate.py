# --- Core Imports ---
import KratosMultiphysics
from KratosMultiphysics.KratosUnittest import TestCase, main

# --- HDF5 Imports ---
import KratosMultiphysics.HDF5Application as HDF5Application


class TestPipedModelPredicates(TestCase):

    @staticmethod
    def __GetModel() -> "tuple[KratosMultiphysics.Model, KratosMultiphysics.ModelPart]":
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("test")
        return model, model_part


    def test_TimeIntervalPredicate(self) -> None:
        model, model_part = self.__GetModel()
        parameters = KratosMultiphysics.Parameters("""[
            {"model_part_name" : "test"}, // <== model part from model
            {},                           // <== process info from model part
            {},                           // <== TIME from process info
            {"interval" : [1.0, 2.0]}     // <== IntervalPredicate == IntervalUtility
        ]""")
        predicate = HDF5Application.TimeIntervalPredicate(parameters)

        for time in (-1.5, 0.0, 0.5, 1.2, 1.8, 2.5, 3.0):
            model_part.ProcessInfo[KratosMultiphysics.TIME] = time
            self.assertEqual(predicate(model), 1.0 <= time and time <= 2.0)


    def test_StepIntervalPredicate(self) -> None:
        model, model_part = self.__GetModel()
        parameters = KratosMultiphysics.Parameters("""[
            {"model_part_name" : "test"}, // <== model part from model
            {},                           // <== process info from model part
            {},                           // <== TIME from process info
            {"interval" : [5, 10]}        // <== DiscreteIntervalPredicate == DiscreteIntervalUtility
        ]""")
        predicate = HDF5Application.StepIntervalPredicate(parameters)

        for step in range(-5, 15):
            model_part.ProcessInfo[KratosMultiphysics.STEP] = step
            self.assertEqual(predicate(model), 5 <= step and step <= 10)


    def test_PeriodicTimeIntervalPredicate(self) -> None:
        model, model_part = self.__GetModel()
        parameters = KratosMultiphysics.Parameters("""[
            {"model_part_name" : "test"}, // <== model part from model
            {},                           // <== process info from model part
            {},                           // <== TIME from process info
            {"mod" : 3.0},                // <== modulo
            {"interval" : [1.0, 2.0]}     // <== IntervalPredicate == IntervalUtility
        ]""")
        predicate = HDF5Application.PeriodicTimeIntervalPredicate(parameters)

        for time, result in ((-4.5, False),
                             (-3.5, False),
                             (-1.5, False),
                             (-0.5, False),
                             (0.0, False),
                             (0.5, False),
                             (1.5, True),
                             (2.5, False),
                             (4.5, True)):
            model_part.ProcessInfo[KratosMultiphysics.TIME] = time
            self.assertEqual(predicate(model), result, msg = str(time))


    def test_PeriodicStepIntervalPredicate(self) -> None:
        model, model_part = self.__GetModel()
        parameters = KratosMultiphysics.Parameters("""[
            {"model_part_name" : "test"}, // <== model part from model
            {},                           // <== process info from model part
            {},                           // <== TIME from process info
            {"mod" : 10},                 // <== modulo
            {"interval" : [5, 10]}        // <== DiscreteIntervalPredicate == DiscreteIntervalUtility
        ]""")
        predicate = HDF5Application.PeriodicStepIntervalPredicate(parameters)

        for step, result in ((-27, False),
                             (-22, False),
                             (-17, False),
                             (-12, False),
                             (-7, False),
                             (-2, False),
                             (0, False),
                             (2, False),
                             (7, True),
                             (12, False),
                             (17, True),
                             (22, False),
                             (27, True)):
            model_part.ProcessInfo[KratosMultiphysics.STEP] = step
            self.assertEqual(predicate(model), result, msg = str(step))


if __name__ == "__main__":
    main()
