import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.DamApplication.impose_water_loads_condition_process import (
    ImposeWaterLoadsConditionProcess,
)


class TestImposeWaterLoadsConditionProcess(KratosUnittest.TestCase):

    def test_constant_hydrostatic_pressure_is_clamped_and_reapplied(self):
        model, model_part = self._CreateModelPart("HydroSurfacePressure3D_Test")
        settings = self._CreateSettings(
            gravity_direction="X",
            specific_weight=10000.0,
            water_level=3.0,
        )

        process = ImposeWaterLoadsConditionProcess(model, settings)

        process.ExecuteBeforeSolutionLoop()
        # pressure = specific_weight * (water_level - x)
        self._CheckPressure(model_part, [30000.0, 10000.0, 0.0])
        self._CheckPressureIsNotFixed(model_part)

        # The legacy Dam process re-applied the constant pressure at every
        # solution step. Verify that the wrapper keeps that behaviour.
        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE, 123.0)
        process.ExecuteInitializeSolutionStep()
        self._CheckPressure(model_part, [30000.0, 10000.0, 0.0])

    def test_table_hydrostatic_pressure_uses_initial_level_then_table(self):
        model, model_part = self._CreateModelPart("HydroSurfacePressure3D_Test")
        model_part.ProcessInfo[KratosDam.TIME_UNIT_CONVERTER] = 2.0

        table = KratosMultiphysics.PiecewiseLinearTable()
        table.AddRow(0.0, 2.0)
        table.AddRow(2.0, 6.0)
        model_part.AddTable(1, table)

        settings = self._CreateSettings(
            gravity_direction="Y",
            specific_weight=10000.0,
            water_level=1.0,
            water_table=1,
        )

        process = ImposeWaterLoadsConditionProcess(model, settings)

        # Before the solution loop the Dam process used Water_level, not the
        # table value at TIME = 0.
        process.ExecuteBeforeSolutionLoop()
        self._CheckPressure(model_part, [10000.0, 0.0, 0.0])

        # TIME = 2.0 s with TIME_UNIT_CONVERTER = 2.0 -> table abscissa 1.0
        # -> water level 4.0.
        model_part.ProcessInfo[KratosMultiphysics.TIME] = 2.0
        process.ExecuteInitializeSolutionStep()
        self._CheckPressure(model_part, [40000.0, 20000.0, 0.0])

    def test_gravity_direction_defaults_to_y(self):
        model, model_part = self._CreateModelPart("HydroSurfacePressure3D_Test")
        settings = self._CreateSettings(
            gravity_direction=None,
            specific_weight=10000.0,
            water_level=3.0,
        )

        process = ImposeWaterLoadsConditionProcess(model, settings)
        process.ExecuteBeforeSolutionLoop()

        self._CheckPressure(model_part, [30000.0, 10000.0, 0.0])

    @staticmethod
    def _CreateModelPart(name):
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart(name)
        model_part.AddNodalSolutionStepVariable(
            KratosMultiphysics.POSITIVE_FACE_PRESSURE
        )
        model_part.ProcessInfo[KratosDam.TIME_UNIT_CONVERTER] = 1.0

        # The coordinates are used by the three gravity directions.
        for node_id, coordinates in enumerate(
            ((0.0, 0.0, 0.0), (2.0, 2.0, 0.0), (5.0, 5.0, 0.0)),
            start=1,
        ):
            model_part.CreateNewNode(node_id, *coordinates)

        return model, model_part

    @staticmethod
    def _CreateSettings(
        gravity_direction,
        specific_weight,
        water_level,
        water_table=0,
    ):
        gravity_direction_entry = ""
        if gravity_direction is not None:
            gravity_direction_entry = (
                '"Gravity_Direction" : "%s",' % gravity_direction
            )

        return KratosMultiphysics.Parameters(
            """
            {
                "model_part_name"                                  : "HydroSurfacePressure3D_Test",
                "variable_name"                                    : "POSITIVE_FACE_PRESSURE",
                "interval"                                         : [0.0, 0.0],
                "Modify"                                           : true,
                %s
                "Reservoir_Bottom_Coordinate_in_Gravity_Direction" : 0.0,
                "Spe_weight"                                       : %.15g,
                "Water_level"                                      : %.15g,
                "Water_Table"                                      : %d
            }
            """
            % (
                gravity_direction_entry,
                specific_weight,
                water_level,
                water_table,
            )
        )

    def _CheckPressure(self, model_part, expected_values):
        for node, expected_value in zip(model_part.Nodes, expected_values):
            self.assertAlmostEqual(
                node.GetSolutionStepValue(
                    KratosMultiphysics.POSITIVE_FACE_PRESSURE
                ),
                expected_value,
            )

    def _CheckPressureIsNotFixed(self, model_part):
        for node in model_part.Nodes:
            self.assertFalse(
                node.IsFixed(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
            )


if __name__ == "__main__":
    KratosUnittest.main()
