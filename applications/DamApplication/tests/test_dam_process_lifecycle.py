import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestDamProcessLifetime(KratosUnittest.TestCase):
    """Assignment-lifetime contracts for Dam scalar-assignment processes.

    These processes were historically implemented with
    ``ApplyConstantScalarValueProcess``, which assigned the value once at
    initialization. After the migration to ``AssignScalarVariableProcess`` the
    wrappers must forward an initial-only interval ``[0.0, 0.0]`` so that the
    value is an initial condition and is not re-applied at later solution steps.
    """

    def _make_model(self, variables, target_submodel_name="target"):
        model = KratosMultiphysics.Model()
        mp = model.CreateModelPart("main")
        mp.ProcessInfo[KratosDam.TIME_UNIT_CONVERTER] = 1.0
        for var in variables:
            mp.AddNodalSolutionStepVariable(var)
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)  # in target model part
        mp.CreateNewNode(2, 5.0, 5.0, 0.0)  # outside target model part
        sub = mp.CreateSubModelPart(target_submodel_name)
        sub.AddNodes([1])
        return model, mp

    def _make_wrapper_settings(self, model_part_name, extra):
        # build {"python_module","Parameters": {...}} as the process_factory would
        params = KratosMultiphysics.Parameters(extra)
        params["model_part_name"].SetString(model_part_name)
        full = KratosMultiphysics.Parameters("""{
            "python_module": "placeholder",
            "kratos_module": "KratosMultiphysics.DamApplication",
            "process_name": "placeholder",
            "Parameters": {}
        }""")
        full["Parameters"] = params
        return full

    def _assert_node_value(
        self, model_part, node_id, variable, expected, tol=1.0e-9, msg=""
    ):
        self.assertAlmostEqual(
            model_part.GetNode(node_id).GetSolutionStepValue(variable),
            expected,
            delta=tol,
            msg=msg,
        )

    def test_uniform_initial_temperature_is_applied_once(self):
        """The uniform initial temperature is an initial condition, not a
        persistent boundary condition: a later solution step must not restore
        the initial value over an evolved temperature field."""
        from KratosMultiphysics.DamApplication.impose_uniform_temperature_process import (
            Factory as UniformTemperatureFactory,
        )

        model, mp = self._make_model(
            [KratosMultiphysics.TEMPERATURE],
            target_submodel_name="INITIALTEMPERATURE_target",
        )
        settings = self._make_wrapper_settings(
            "main.INITIALTEMPERATURE_target",
            r"""
            {
                "model_part_name": "",
                "variable_name": "TEMPERATURE",
                "interval": [0.0, 0.0],
                "constrained": false,
                "value": 12.712,
                "table": 0
            }""",
        )
        process = UniformTemperatureFactory(settings, model)

        process.ExecuteBeforeSolutionLoop()
        self._assert_node_value(mp, 1, KratosMultiphysics.TEMPERATURE, 12.712)

        # emulate thermal evolution (solved field)
        mp.GetNode(1).SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 8.0)
        mp.ProcessInfo[KratosMultiphysics.TIME] = 2.0
        mp.CloneTimeStep(2.0)
        process.ExecuteInitializeSolutionStep()

        # the later step must NOT overwrite the evolved temperature
        self._assert_node_value(
            mp,
            1,
            KratosMultiphysics.TEMPERATURE,
            8.0,
            msg="initial temperature was re-applied in a later step",
        )

    def test_thermal_parameters_are_initialized_once(self):
        """DENSITY / CONDUCTIVITY / SPECIFIC_HEAT are applied once at
        initialization; a later step must not restore the initialized value over
        a value updated by another process/model."""
        from KratosMultiphysics.DamApplication.impose_thermal_parameters_scalar_value_process import (
            Factory as ThermalParametersFactory,
        )

        model, mp = self._make_model(
            [
                KratosMultiphysics.CONDUCTIVITY,
                KratosMultiphysics.SPECIFIC_HEAT,
                KratosMultiphysics.DENSITY,
            ],
            target_submodel_name="ThermalParameters_target",
        )
        settings = self._make_wrapper_settings(
            "main.ThermalParameters_target",
            r"""
            {
                "model_part_name": "",
                "variable_name": "THERMAL_PARAMETERS",
                "interval": [0.0, 0.0],
                "ThermalDensity": 2400,
                "Conductivity": 2.2,
                "SpecificHeat": 1000.0
            }""",
        )
        process = ThermalParametersFactory(settings, model)

        process.ExecuteBeforeSolutionLoop()
        self._assert_node_value(mp, 1, KratosMultiphysics.CONDUCTIVITY, 2.2)
        self._assert_node_value(mp, 1, KratosMultiphysics.DENSITY, 2400.0)
        self._assert_node_value(mp, 1, KratosMultiphysics.SPECIFIC_HEAT, 1000.0)

        # emulate a later update of the conductivity by another process/model
        mp.GetNode(1).SetSolutionStepValue(KratosMultiphysics.CONDUCTIVITY, 5.0)
        mp.ProcessInfo[KratosMultiphysics.TIME] = 2.0
        mp.CloneTimeStep(2.0)
        process.ExecuteInitializeSolutionStep()

        # the later step must NOT overwrite the updated conductivity
        self._assert_node_value(
            mp,
            1,
            KratosMultiphysics.CONDUCTIVITY,
            5.0,
            msg="thermal parameter was re-applied in a later step",
        )

    def test_uniform_face_heat_flux_is_not_reapplied(self):
        """The uniform face heat flux is assigned once at initialization and
        remains unchanged unless another process/model modifies it."""
        from KratosMultiphysics.DamApplication.impose_face_heat_flux_process import (
            Factory as FaceHeatFluxFactory,
        )

        model, mp = self._make_model(
            [KratosMultiphysics.FACE_HEAT_FLUX],
            target_submodel_name="UniformFlux3D_target",
        )
        settings = self._make_wrapper_settings(
            "main.UniformFlux3D_target",
            r"""
            {
                "model_part_name": "",
                "variable_name": "FACE_HEAT_FLUX",
                "interval": [0.0, 0.0],
                "constrained": false,
                "value": 10.0,
                "table": 0
            }""",
        )
        process = FaceHeatFluxFactory(settings, model)

        process.ExecuteBeforeSolutionLoop()
        self._assert_node_value(mp, 1, KratosMultiphysics.FACE_HEAT_FLUX, 10.0)

        # emulate a later update of the flux
        mp.GetNode(1).SetSolutionStepValue(KratosMultiphysics.FACE_HEAT_FLUX, 20.0)
        mp.ProcessInfo[KratosMultiphysics.TIME] = 2.0
        mp.CloneTimeStep(2.0)
        process.ExecuteInitializeSolutionStep()

        # the later step must NOT restore the initialized flux
        self._assert_node_value(
            mp,
            1,
            KratosMultiphysics.FACE_HEAT_FLUX,
            20.0,
            msg="uniform face heat flux was re-applied in a later step",
        )


if __name__ == "__main__":
    KratosUnittest.main()
