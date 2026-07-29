import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.DamApplication.apply_load_vector_dam_table_process import (
    ApplyLoadVectorDamTableProcess,
)
from KratosMultiphysics.DamApplication.apply_load_vector_dam_process import (
    ApplyLoadVectorDamProcess,
)


class TestApplyLoadVectorDamProcesses(KratosUnittest.TestCase):

    def test_load_does_not_fix_components(self):
        model, model_part = self._CreateModelPart()
        process = ApplyLoadVectorDamProcess(model, self._CreateSettings())

        process.ExecuteBeforeSolutionLoop()

        self._CheckNodalValues(model_part, [3.0, -6.0, 1.5])
        self._CheckComponentsAreFree(model_part)

    def test_constant_table_load_does_not_fix_components(self):
        model, model_part = self._CreateModelPart()
        settings = self._CreateSettings()
        settings.AddEmptyValue("table").SetInt(0)
        process = ApplyLoadVectorDamTableProcess(
            model,
            settings,
        )

        process.ExecuteBeforeSolutionLoop()

        self._CheckNodalValues(model_part, [3.0, -6.0, 1.5])
        self._CheckComponentsAreFree(model_part)

    def test_table_load_does_not_fix_components(self):
        model, model_part = self._CreateModelPart()
        table = KratosMultiphysics.PiecewiseLinearTable()
        table.AddRow(0.0, 2.0)
        table.AddRow(1.0, 4.0)
        model_part.AddTable(1, table)

        settings = self._CreateSettings()
        settings.AddEmptyValue("table").SetInt(1)
        process = ApplyLoadVectorDamTableProcess(
            model,
            settings,
        )

        process.ExecuteBeforeSolutionLoop()
        self._CheckComponentsAreFree(model_part)

        model_part.CloneTimeStep(1.0)
        process.ExecuteInitializeSolutionStep()

        self._CheckNodalValues(model_part, [4.0, 4.0, 4.0])
        self._CheckComponentsAreFree(model_part)

    @staticmethod
    def _CreateModelPart():
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.ProcessInfo[KratosDam.TIME_UNIT_CONVERTER] = 1.0

        for node_id, coordinates in enumerate(
            ((0.0, 0.0, 0.0), (1.0, 2.0, 3.0)),
            start=1,
        ):
            node = model_part.CreateNewNode(node_id, *coordinates)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z)

        return model, model_part

    @staticmethod
    def _CreateSettings():
        return KratosMultiphysics.Parameters(
            """
            {
                "model_part_name" : "main",
                "variable_name"   : "DISPLACEMENT",
                "modulus"         : 3.0,
                "direction"       : [1.0, -2.0, 0.5]
            }
            """
        )

    def _CheckNodalValues(self, model_part, expected_values):
        variables = (
            KratosMultiphysics.DISPLACEMENT_X,
            KratosMultiphysics.DISPLACEMENT_Y,
            KratosMultiphysics.DISPLACEMENT_Z,
        )
        for node in model_part.Nodes:
            for variable, expected_value in zip(variables, expected_values):
                self.assertAlmostEqual(
                    node.GetSolutionStepValue(variable),
                    expected_value,
                )

    def _CheckComponentsAreFree(self, model_part):
        variables = (
            KratosMultiphysics.DISPLACEMENT_X,
            KratosMultiphysics.DISPLACEMENT_Y,
            KratosMultiphysics.DISPLACEMENT_Z,
        )
        for node in model_part.Nodes:
            for variable in variables:
                self.assertFalse(node.IsFixed(variable))


if __name__ == "__main__":
    KratosUnittest.main()
