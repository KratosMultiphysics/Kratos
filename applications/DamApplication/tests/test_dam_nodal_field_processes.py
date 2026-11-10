import json
import os
import tempfile

import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.DamApplication.impose_nodal_reference_temperature_process import (
    Factory as ReferenceTemperatureFactory,
)
from KratosMultiphysics.DamApplication.impose_input_table_nodal_young_modulus_process import (
    Factory as InputTableYoungModulusFactory,
)


class TestDamNodalFieldProcesses(KratosUnittest.TestCase):
    """End-to-end checks of the Dam nodal-field producers through their Python
    wrappers (the path used by the GiD interface and ProjectParameters)."""

    def setUp(self):
        self.model = KratosMultiphysics.Model()
        self.model_part = self.model.CreateModelPart("main")
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.YOUNG_MODULUS)
        self.model_part.AddNodalSolutionStepVariable(KratosDam.NODAL_YOUNG_MODULUS)
        coordinates = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
        for node_id, node_coordinates in enumerate(coordinates, start=1):
            self.model_part.CreateNewNode(node_id, *node_coordinates)
        self.model_part.CreateNewProperties(1)
        self.temporary_directory = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temporary_directory.cleanup()

    def _WriteTableFile(self, rows):
        file_path = os.path.join(self.temporary_directory.name, "nodal_values.txt")
        with open(file_path, "w") as table_file:
            for row in rows:
                table_file.write("{} {}\n".format(row[0], row[1]))
        return file_path

    def _ReferenceTemperatureSettings(self, input_file_name, initial_value=11.0):
        return KratosMultiphysics.Parameters(
            """
            {
                "Parameters" : {
                    "model_part_name" : "main",
                    "variable_name"   : "TEMPERATURE",
                    "initial_value"   : %.15g,
                    "input_file_name" : %s
                }
            }
            """
            % (initial_value, json.dumps(input_file_name))
        )

    def _InputTableYoungSettings(self, input_file_name, initial_value, variable_name):
        return KratosMultiphysics.Parameters(
            """
            {
                "Parameters" : {
                    "model_part_name" : "main",
                    "variable_name"   : %s,
                    "initial_value"   : %.15g,
                    "min_value"       : 0.0,
                    "max_value"       : 0.0,
                    "input_file_name" : %s
                }
            }
            """
            % (json.dumps(variable_name), initial_value, json.dumps(input_file_name))
        )

    def _NodalValues(self, variable):
        return [node.GetSolutionStepValue(variable) for node in self.model_part.Nodes]

    def test_reference_temperature_initial_value_without_input_file(self):
        process = ReferenceTemperatureFactory(
            self._ReferenceTemperatureSettings(""), self.model
        )

        process.ExecuteBeforeSolutionLoop()

        self.assertEqual(self._NodalValues(KratosMultiphysics.TEMPERATURE), [11.0] * 3)

        # The initial value is re-applied at every solution step.
        self.model_part.CloneTimeStep(1.0)
        process.ExecuteInitializeSolutionStep()
        self.assertEqual(self._NodalValues(KratosMultiphysics.TEMPERATURE), [11.0] * 3)

    def test_reference_temperature_table_values_by_node_id(self):
        file_path = self._WriteTableFile(((1, 100.0), (2, 250.0), (3, -30.0)))
        process = ReferenceTemperatureFactory(
            self._ReferenceTemperatureSettings(file_path), self.model
        )

        process.ExecuteBeforeSolutionLoop()
        self.assertEqual(
            self._NodalValues(KratosMultiphysics.TEMPERATURE), [100.0, 250.0, -30.0]
        )

        self.model_part.CloneTimeStep(1.0)
        process.ExecuteInitializeSolutionStep()
        self.assertEqual(
            self._NodalValues(KratosMultiphysics.TEMPERATURE), [100.0, 250.0, -30.0]
        )

    def test_input_table_nodal_young_modulus_initial_value_installs_accessor(self):
        process = InputTableYoungModulusFactory(
            self._InputTableYoungSettings("- No file", 3.0e7, "NODAL_YOUNG_MODULUS"),
            self.model,
        )

        process.ExecuteBeforeSolutionLoop()

        self.assertEqual(
            self._NodalValues(KratosMultiphysics.YOUNG_MODULUS), [3.0e7] * 3
        )
        # The legacy alias is never written.
        self.assertEqual(self._NodalValues(KratosDam.NODAL_YOUNG_MODULUS), [0.0] * 3)
        self.assertTrue(
            self.model_part.GetProperties()[1].HasAccessor(
                KratosMultiphysics.YOUNG_MODULUS
            )
        )

        # Repeated initialization must be idempotent.
        process.ExecuteBeforeSolutionLoop()
        self.assertEqual(
            self._NodalValues(KratosMultiphysics.YOUNG_MODULUS), [3.0e7] * 3
        )
        self.assertTrue(
            self.model_part.GetProperties()[1].HasAccessor(
                KratosMultiphysics.YOUNG_MODULUS
            )
        )

    def test_input_table_nodal_young_modulus_table_values_by_node_id(self):
        file_path = self._WriteTableFile(((1, 1.0e7), (2, 2.0e7), (3, 4.0e7)))
        process = InputTableYoungModulusFactory(
            self._InputTableYoungSettings(file_path, 0.0, "YOUNG_MODULUS"), self.model
        )

        process.ExecuteBeforeSolutionLoop()
        process.ExecuteInitializeSolutionStep()

        self.assertEqual(
            self._NodalValues(KratosMultiphysics.YOUNG_MODULUS),
            [1.0e7, 2.0e7, 4.0e7],
        )
        self.assertEqual(self._NodalValues(KratosDam.NODAL_YOUNG_MODULUS), [0.0] * 3)
        self.assertTrue(
            self.model_part.GetProperties()[1].HasAccessor(
                KratosMultiphysics.YOUNG_MODULUS
            )
        )


if __name__ == "__main__":
    KratosUnittest.main()
