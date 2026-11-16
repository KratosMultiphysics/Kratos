import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.PoromechanicsApplication as KratosPoro

## This process sets the value of water loads.

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeWaterLoadsConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ImposeWaterLoadsConditionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]
        self.model_part = model_part

        self.components_process_list = []
        self.hydrostatic_processes = []
        self.uplift_coupling_data = []
        # The Poromechanics constant hydrostatic process only sets the pressure
        # during ExecuteInitialize. The legacy Dam process re-applied the same
        # constant pressure at every solution step, so keep doing that here.
        self.constant_hydrostatic_processes = []

        if ("HydroLinePressure2D" in settings["model_part_name"].GetString()) or ("HydroSurfacePressure3D" in settings["model_part_name"].GetString()):

            hydrostatic_settings = self._CreateHydrostaticSettings(settings)

            if hydrostatic_settings["table"].GetInt() == 0:
                process = KratosPoro.ApplyConstantHydrostaticPressureProcess(model_part, hydrostatic_settings)
                self.constant_hydrostatic_processes.append(process)
            else:
                process = KratosPoro.ApplyHydrostaticPressureTableProcess(model_part, hydrostatic_settings)

            self.hydrostatic_processes.append(process)
            self.components_process_list.append(process)

        if ("StraightUpliftLinePressure2D" in settings["model_part_name"].GetString()) or ("StraightUpliftSurfacePressure3D" in settings["model_part_name"].GetString()):

            joint_model_part = Model["MainModelPart.Parts_" + settings["joint_group_name"].GetString()]
            uplift_process = KratosDam.DamUpliftConditionLoadProcess(
                model_part, joint_model_part, self._UpliftProcessSettings(settings))
            self.components_process_list.append(uplift_process)

            coupling_settings = self._CreateUpliftCouplingSettings(settings)
            if coupling_settings["enabled"].GetBool():
                self.uplift_coupling_data.append(
                    {
                        "process": uplift_process,
                        "model_part": model_part,
                        "variable": KratosMultiphysics.KratosGlobals.GetVariable(
                            settings["variable_name"].GetString()),
                        "settings": coupling_settings,
                    }
                )

        if "CircularUpliftSurfacePressure3D" in settings["model_part_name"].GetString():

            joint_model_part = Model["MainModelPart.Parts_" + settings["joint_group_name"].GetString()]
            self.components_process_list.append(KratosDam.DamUpliftCircularConditionLoadProcess(model_part, joint_model_part, settings))

        if ("HydroDynamicLinePressure2D" in settings["model_part_name"].GetString()) or ("HydroDynamicSurfacePressure3D" in settings["model_part_name"].GetString()):

            self.components_process_list.append(KratosDam.DamWestergaardConditionLoadProcess(model_part, settings))

    @staticmethod
    def _UpliftProcessSettings(settings):
        # The hydro-mechanical coupling options are consumed by the Python
        # coupling controller and are not part of the C++ uplift process
        # schema. Pass a copy without them.
        if not settings.Has("hydro_mechanical_coupling"):
            return settings
        uplift_settings = KratosMultiphysics.Parameters(settings.WriteJsonString())
        uplift_settings.RemoveValue("hydro_mechanical_coupling")
        return uplift_settings

    @staticmethod
    def _CreateUpliftCouplingSettings(settings):
        """Return the coupling settings with defaults applied.

        The coupling is enabled by default. It requires a quasi-static
        Newton-Raphson mechanical solve (dynamic and arc-length are rejected)
        and all straight uplift processes must share identical settings.
        """
        coupling_settings = KratosMultiphysics.Parameters(
            """
            {
                "enabled"                         : true,
                "max_iterations"                  : 50,
                "pressure_relative_tolerance"     : 1.0e-8,
                "pressure_absolute_tolerance"     : 1.0e-4,
                "joint_width_absolute_tolerance"  : 1.0e-7,
                "front_relative_tolerance"        : 1.0e-6,
                "relaxation_factor"               : 1.0,
                "echo_level"                      : 0
            }
            """)
        if settings.Has("hydro_mechanical_coupling"):
            user_settings = settings["hydro_mechanical_coupling"]
            user_settings.ValidateAndAssignDefaults(coupling_settings)
            coupling_settings = user_settings
        return coupling_settings

    def GetUpliftCouplingData(self):
        """Return the description of the straight-uplift processes that must be
        re-evaluated during the hydro-mechanical coupling iterations."""
        return list(self.uplift_coupling_data)

    @staticmethod
    def _CreateHydrostaticSettings(settings):
        # Adapt the Dam input schema to the Poromechanics process without
        # changing the Dam project files. The legacy Dam process read the
        # reservoir bottom coordinate but used Water_level directly, so
        # Water_level is the Poromechanics reference_coordinate.
        hydrostatic_settings = KratosMultiphysics.Parameters("{}")
        hydrostatic_settings.AddValue("model_part_name", settings["model_part_name"])
        hydrostatic_settings.AddValue("variable_name", settings["variable_name"])
        hydrostatic_settings.AddEmptyValue("is_fixed").SetBool(False)

        gravity_direction = settings["Gravity_Direction"].GetString() if settings.Has("Gravity_Direction") else "Y"
        if gravity_direction == "X":
            direction = 0
        elif gravity_direction == "Y":
            direction = 1
        else:
            direction = 2
        hydrostatic_settings.AddEmptyValue("gravity_direction").SetInt(direction)

        water_level = settings["Water_level"].GetDouble() if settings.Has("Water_level") else 0.0
        specific_weight = settings["Spe_weight"].GetDouble() if settings.Has("Spe_weight") else 0.0
        water_table = settings["Water_Table"].GetInt() if settings.Has("Water_Table") else 0

        hydrostatic_settings.AddEmptyValue("reference_coordinate").SetDouble(water_level)
        hydrostatic_settings.AddEmptyValue("specific_weight").SetDouble(specific_weight)
        hydrostatic_settings.AddEmptyValue("table").SetInt(water_table)

        return hydrostatic_settings

    def ExecuteBeforeSolutionLoop(self):

        # The Poromechanics processes initialize hydrostatic pressures in
        # ExecuteInitialize; the legacy Dam process did it before the loop.
        for process in self.hydrostatic_processes:
            process.ExecuteInitialize()

        for component in self.components_process_list:
            component.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):

        for process in self.constant_hydrostatic_processes:
            process.ExecuteInitialize()

        for process in self.components_process_list:
            process.ExecuteInitializeSolutionStep()
