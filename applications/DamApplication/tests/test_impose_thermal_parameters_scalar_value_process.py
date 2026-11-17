import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.process_factory import KratosProcessFactory

from KratosMultiphysics.DamApplication.impose_thermal_parameters_scalar_value_process import (
    ImposeThemalParametersScalarValueProcess,
)


class TestImposeThermalParametersScalarValueProcess(KratosUnittest.TestCase):
    """Focused checks of the DENSITY / CONDUCTIVITY / SPECIFIC_HEAT assignment
    wrapper, including its initialization-only behavior."""

    sentinel_values = (111.0, 222.0, 333.0)

    def setUp(self):
        self.model = KratosMultiphysics.Model()
        self.model_part = self.model.CreateModelPart("main")
        for variable in self._Variables():
            self.model_part.AddNodalSolutionStepVariable(variable)
        for node_id in (1, 2):
            self.model_part.CreateNewNode(node_id, float(node_id), 0.0, 0.0)

    @staticmethod
    def _Variables():
        return (
            KratosMultiphysics.DENSITY,
            KratosMultiphysics.CONDUCTIVITY,
            KratosMultiphysics.SPECIFIC_HEAT,
        )

    @staticmethod
    def _CreateSettings(
        thermal_density=2400.0,
        conductivity=2.2,
        specific_heat=1000.0,
        interval=(0.0, 0.0),
    ):
        return KratosMultiphysics.Parameters(
            """
            {
                "model_part_name" : "main",
                "variable_name"   : "THERMAL_PARAMETERS",
                "interval"        : [%.15g, %.15g],
                "ThermalDensity"  : %.15g,
                "Conductivity"    : %.15g,
                "SpecificHeat"    : %.15g
            }
            """
            % (
                interval[0],
                interval[1],
                thermal_density,
                conductivity,
                specific_heat,
            )
        )

    def _SetSentinels(self):
        for node in self.model_part.Nodes:
            for variable, value in zip(self._Variables(), self.sentinel_values):
                node.SetSolutionStepValue(variable, value)

    def _NodalValues(self, variable):
        return [node.GetSolutionStepValue(variable) for node in self.model_part.Nodes]

    def _CheckAllNodes(self, variable, expected):
        for value in self._NodalValues(variable):
            self.assertAlmostEqual(value, expected)

    def test_factory_assigns_all_thermal_parameters(self):
        process = ImposeThemalParametersScalarValueProcess(
            self.model,
            self._CreateSettings(thermal_density=2400.0, conductivity=2.2, specific_heat=1000.0),
        )

        process.ExecuteBeforeSolutionLoop()

        self._CheckAllNodes(KratosMultiphysics.DENSITY, 2400.0)
        self._CheckAllNodes(KratosMultiphysics.CONDUCTIVITY, 2.2)
        self._CheckAllNodes(KratosMultiphysics.SPECIFIC_HEAT, 1000.0)

    def test_parameter_order_is_preserved(self):
        process = ImposeThemalParametersScalarValueProcess(
            self.model, self._CreateSettings()
        )

        self.assertEqual(
            [p.variable.Name() for p in process.components_process_list],
            ["DENSITY", "CONDUCTIVITY", "SPECIFIC_HEAT"],
        )

    def test_zero_values_are_skipped(self):
        self._SetSentinels()
        process = ImposeThemalParametersScalarValueProcess(
            self.model,
            self._CreateSettings(
                thermal_density=0.0, conductivity=2.2, specific_heat=0.0
            ),
        )

        self.assertEqual(len(process.components_process_list), 1)

        process.ExecuteBeforeSolutionLoop()

        # Only the nonzero parameter is assigned.
        self._CheckAllNodes(
            KratosMultiphysics.DENSITY, self.sentinel_values[0]
        )
        self._CheckAllNodes(KratosMultiphysics.CONDUCTIVITY, 2.2)
        self._CheckAllNodes(
            KratosMultiphysics.SPECIFIC_HEAT, self.sentinel_values[2]
        )

    def test_all_zero_values_create_no_processes(self):
        self._SetSentinels()
        process = ImposeThemalParametersScalarValueProcess(
            self.model,
            self._CreateSettings(
                thermal_density=0.0, conductivity=0.0, specific_heat=0.0
            ),
        )

        self.assertEqual(len(process.components_process_list), 0)

        process.ExecuteBeforeSolutionLoop()

        for variable, sentinel in zip(self._Variables(), self.sentinel_values):
            self._CheckAllNodes(variable, sentinel)

    def test_negative_values_are_assigned(self):
        process = ImposeThemalParametersScalarValueProcess(
            self.model,
            self._CreateSettings(
                thermal_density=-2400.0, conductivity=-2.2, specific_heat=-1000.0
            ),
        )

        self.assertEqual(len(process.components_process_list), 3)

        process.ExecuteBeforeSolutionLoop()

        self._CheckAllNodes(KratosMultiphysics.DENSITY, -2400.0)
        self._CheckAllNodes(KratosMultiphysics.CONDUCTIVITY, -2.2)
        self._CheckAllNodes(KratosMultiphysics.SPECIFIC_HEAT, -1000.0)

    def test_values_are_applied_once_at_time_zero(self):
        process = ImposeThemalParametersScalarValueProcess(
            self.model, self._CreateSettings()
        )

        process.ExecuteBeforeSolutionLoop()
        self._CheckAllNodes(KratosMultiphysics.CONDUCTIVITY, 2.2)

        # Emulate a later update of the parameters by another process/model.
        for variable in self._Variables():
            for node in self.model_part.Nodes:
                node.SetSolutionStepValue(variable, 7.0)
        self.model_part.ProcessInfo[KratosMultiphysics.TIME] = 2.0
        self.model_part.CloneTimeStep(2.0)
        process.ExecuteInitializeSolutionStep()

        # The interval [0.0, 0.0] keeps the assignment at the initial step only.
        for variable in self._Variables():
            self._CheckAllNodes(variable, 7.0)

    def test_nothing_is_assigned_when_starting_after_time_zero(self):
        self._SetSentinels()
        self.model_part.ProcessInfo[KratosMultiphysics.TIME] = 5.0
        process = ImposeThemalParametersScalarValueProcess(
            self.model, self._CreateSettings()
        )

        process.ExecuteBeforeSolutionLoop()

        # Documented legacy behavior: the fixed interval [0.0, 0.0] is not
        # active when the simulation starts after TIME = 0.
        for variable, sentinel in zip(self._Variables(), self.sentinel_values):
            self._CheckAllNodes(variable, sentinel)

    def test_incoming_interval_parameter_is_ignored(self):
        # The ProjectParameters interval is not forwarded; the component
        # processes always use [0.0, 0.0].
        process = ImposeThemalParametersScalarValueProcess(
            self.model, self._CreateSettings(interval=(5.0, 6.0))
        )

        process.ExecuteBeforeSolutionLoop()

        self._CheckAllNodes(KratosMultiphysics.DENSITY, 2400.0)
        self._CheckAllNodes(KratosMultiphysics.CONDUCTIVITY, 2.2)
        self._CheckAllNodes(KratosMultiphysics.SPECIFIC_HEAT, 1000.0)

    def test_process_factory_constructs_the_process(self):
        process_list = KratosMultiphysics.Parameters(
            """
            [{
                "python_module" : "impose_thermal_parameters_scalar_value_process",
                "kratos_module" : "KratosMultiphysics.DamApplication",
                "process_name"  : "ImposeThermalParametersScalarValueProcess",
                "Parameters"    : {
                    "model_part_name" : "main",
                    "variable_name"   : "THERMAL_PARAMETERS",
                    "interval"        : [0.0, 0.0],
                    "ThermalDensity"  : 2400.0,
                    "Conductivity"    : 2.2,
                    "SpecificHeat"    : 1000.0
                }
            }]
            """
        )

        processes = KratosProcessFactory(self.model).ConstructListOfProcesses(
            process_list
        )

        self.assertEqual(len(processes), 1)
        processes[0].ExecuteBeforeSolutionLoop()
        self._CheckAllNodes(KratosMultiphysics.DENSITY, 2400.0)
        self._CheckAllNodes(KratosMultiphysics.CONDUCTIVITY, 2.2)
        self._CheckAllNodes(KratosMultiphysics.SPECIFIC_HEAT, 1000.0)


if __name__ == "__main__":
    KratosUnittest.main()
