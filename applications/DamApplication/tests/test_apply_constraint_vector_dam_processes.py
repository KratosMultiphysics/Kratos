import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.DamApplication.apply_constraint_vector_dam_table_process import (
    ApplyConstraintVectorDamTableProcess,
)


class TestApplyConstraintVectorDamProcesses(KratosUnittest.TestCase):

    def test_constant_constraints_set_values_and_fixity(self):
        model, model_part = self._CreateModelPart()
        process = ApplyConstraintVectorDamTableProcess(
            model,
            self._CreateSettings(
                constrained=(True, False, True),
                value=(1.0, 2.0, 3.0),
            ),
        )

        # The three components are always created, in X, Y, Z order.
        self.assertEqual(len(process.components_process_list), 3)
        self.assertEqual(
            [p.variable.Name() for p in process.components_process_list],
            ["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z"],
        )

        process.ExecuteBeforeSolutionLoop()

        self._CheckNodalValues(model_part, [1.0, 2.0, 3.0])
        self._CheckFixity(model_part, [True, False, True])

        # The historical Dam wrapper does not forward ExecuteFinalizeSolutionStep,
        # so the prescribed components stay fixed after the step.
        process.ExecuteInitializeSolutionStep()
        self._CheckNodalValues(model_part, [1.0, 2.0, 3.0])
        self._CheckFixity(model_part, [True, False, True])

    def test_zero_valued_components_are_still_created(self):
        model, model_part = self._CreateModelPart()
        process = ApplyConstraintVectorDamTableProcess(
            model,
            self._CreateSettings(
                constrained=(True, False, False),
                value=(0.0, 0.0, 0.0),
            ),
        )

        # Unlike the load wrappers, constraints are never filtered by value.
        self.assertEqual(len(process.components_process_list), 3)

        process.ExecuteBeforeSolutionLoop()

        self._CheckNodalValues(model_part, [0.0, 0.0, 0.0])
        self._CheckFixity(model_part, [True, False, False])

    def test_table_constraints_use_initial_value_then_table_value(self):
        model, model_part = self._CreateModelPart()
        table_x = KratosMultiphysics.PiecewiseLinearTable()
        table_x.AddRow(0.0, 2.0)
        table_x.AddRow(2.0, 4.0)
        model_part.AddTable(1, table_x)
        table_z = KratosMultiphysics.PiecewiseLinearTable()
        table_z.AddRow(0.0, -1.0)
        table_z.AddRow(2.0, -3.0)
        model_part.AddTable(2, table_z)

        process = ApplyConstraintVectorDamTableProcess(
            model,
            self._CreateSettings(
                constrained=(True, True, False),
                value=(100.0, 5.0, 200.0),
                value_table=(1, 0, 2),
            ),
        )

        # Before the solution loop the prescribed initial values are written.
        process.ExecuteBeforeSolutionLoop()
        self._CheckNodalValues(model_part, [100.0, 5.0, 200.0])
        self._CheckFixity(model_part, [True, True, False])

        # After the first step each table component is evaluated at the current
        # time; the constant component keeps its value.
        model_part.CloneTimeStep(2.0)
        process.ExecuteInitializeSolutionStep()
        self._CheckNodalValues(model_part, [4.0, 5.0, -3.0])
        self._CheckFixity(model_part, [True, True, False])

    def test_table_constraint_respects_time_unit_converter(self):
        model, model_part = self._CreateModelPart()
        model_part.ProcessInfo[KratosDam.TIME_UNIT_CONVERTER] = 2.0
        table = KratosMultiphysics.PiecewiseLinearTable()
        table.AddRow(0.0, 0.0)
        table.AddRow(1.0, 10.0)
        model_part.AddTable(1, table)

        process = ApplyConstraintVectorDamTableProcess(
            model,
            self._CreateSettings(
                constrained=(True, False, False),
                value=(0.0, 0.0, 0.0),
                value_table=(1, 0, 0),
            ),
        )

        process.ExecuteBeforeSolutionLoop()
        # TIME = 2.0 s with TIME_UNIT_CONVERTER = 2.0 -> table abscissa 1.0
        model_part.CloneTimeStep(2.0)
        process.ExecuteInitializeSolutionStep()

        self._CheckNodalValues(model_part, [10.0, 0.0, 0.0])
        self._CheckFixity(model_part, [True, False, False])

    def test_interval_parameter_is_not_propagated_to_component_processes(self):
        model, model_part = self._CreateModelPart()
        settings = self._CreateSettings(
            constrained=(True, True, True),
            value=(1.0, 2.0, 3.0),
        )
        # The Dam wrapper deliberately ignores its own interval. The component
        # process therefore keeps the standard default interval [0.0, 1e30].
        settings.AddValue("interval", KratosMultiphysics.Parameters("[-10.0, -5.0]"))

        process = ApplyConstraintVectorDamTableProcess(model, settings)
        process.ExecuteBeforeSolutionLoop()

        # TIME = 0 lies outside the configured interval, yet the constraint is
        # applied because the interval is not forwarded.
        self._CheckNodalValues(model_part, [1.0, 2.0, 3.0])
        self._CheckFixity(model_part, [True, True, True])

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
    def _CreateSettings(constrained, value, value_table=(0, 0, 0)):
        constrained_flags = ["true" if flag else "false" for flag in constrained]
        return KratosMultiphysics.Parameters(
            """
            {
                "model_part_name" : "main",
                "variable_name"   : "DISPLACEMENT",
                "constrained"     : [%s, %s, %s],
                "value"           : [%.15g, %.15g, %.15g],
                "Value_Table"     : [%d, %d, %d]
            }
            """
            % (
                constrained_flags[0],
                constrained_flags[1],
                constrained_flags[2],
                value[0],
                value[1],
                value[2],
                value_table[0],
                value_table[1],
                value_table[2],
            )
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
                    node.GetSolutionStepValue(variable), expected_value
                )

    def _CheckFixity(self, model_part, expected_fixity):
        variables = (
            KratosMultiphysics.DISPLACEMENT_X,
            KratosMultiphysics.DISPLACEMENT_Y,
            KratosMultiphysics.DISPLACEMENT_Z,
        )
        for node in model_part.Nodes:
            for variable, expected_fixed in zip(variables, expected_fixity):
                self.assertEqual(node.IsFixed(variable), expected_fixed)


if __name__ == "__main__":
    KratosUnittest.main()
