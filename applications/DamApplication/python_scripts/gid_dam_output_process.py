import KratosMultiphysics as KM
from KratosMultiphysics.gid_output_process import GiDOutputProcess


class GiDDamOutputProcess(GiDOutputProcess):
    """Compatibility adapter for the legacy Dam GiD output configuration.

    The legacy Dam settings are translated to the standard GiDOutputProcess
    configuration. Result writing is fully delegated to the standard process.
    """

    legacy_defaults = KM.Parameters('''{
        "result_file_configuration": {
            "gidpost_flags": {
                "GiDPostMode": "GiD_PostBinary",
                "WriteDeformedMeshFlag": "WriteUndeformed",
                "WriteConditionsFlag": "WriteElementsOnly",
                "MultiFileFlag": "SingleFile"
            },
            "output_control_type": "time_s",
            "output_interval": 1.0,
            "start_output_results": 0,
            "body_output": true,
            "node_output": false,
            "skin_output": false,
            "plane_output": [],
            "nodal_results": [],
            "nodal_nonhistorical_results": [],
            "nodal_flags_results": [],
            "gauss_point_results": [],
            "additional_list_files": []
        },
        "point_data_configuration": []
    }''')

    _time_unit_factors = {
        "time_s": 1.0,
        "time_h": 3600.0,
        "time_d": 86400.0,
        "time_w": 604800.0
    }

    def __init__(self, model_part, file_name, start_time, param=None):
        legacy_param = self.legacy_defaults.Clone() if param is None else param.Clone()
        self._prepare_legacy_parameters(legacy_param)
        standard_param = self._translate_parameters(model_part, start_time, legacy_param)

        super().__init__(model_part, file_name, standard_param)

        # GiDOutputProcess formats result labels through this private hook.
        # Overriding it here keeps the legacy output-unit labels without
        # duplicating the result-writing implementation.
        self._GiDOutputProcess__get_pretty_time = self._get_legacy_pretty_time
        self._dam_next_output = 0.0

    def _prepare_legacy_parameters(self, settings):
        if settings.Has("result_file_configuration"):
            result_file_configuration = settings["result_file_configuration"]

            if result_file_configuration.Has("output_frequency"):
                if not result_file_configuration.Has("output_interval"):
                    result_file_configuration.AddEmptyValue("output_interval")
                    result_file_configuration["output_interval"].SetDouble(
                        result_file_configuration["output_frequency"].GetDouble())
                result_file_configuration.RemoveValue("output_frequency")

            if result_file_configuration.Has("write_properties_id"):
                result_file_configuration.RemoveValue("write_properties_id")
                KM.Logger.PrintWarning(
                    "GiDDamOutputProcess",
                    "'write_properties_id' is not supported by the GiD output and is ignored."
                )

        settings.ValidateAndAssignDefaults(self.legacy_defaults)

        if settings.Has("result_file_configuration"):
            self._MoveNodalYoungModulusToNonhistorical(settings["result_file_configuration"])

    @staticmethod
    def _MoveNodalYoungModulusToNonhistorical(result_file_configuration):
        variable_name = "NODAL_YOUNG_MODULUS"

        if not result_file_configuration.Has("nodal_results"):
            return

        nodal_results = result_file_configuration["nodal_results"].GetStringArray()
        if variable_name not in nodal_results:
            return

        result_file_configuration["nodal_results"].SetStringArray(
            [name for name in nodal_results if name != variable_name])

        if not result_file_configuration.Has("nodal_nonhistorical_results"):
            result_file_configuration.AddEmptyArray("nodal_nonhistorical_results")
        nonhistorical_results = result_file_configuration["nodal_nonhistorical_results"].GetStringArray()
        if variable_name not in nonhistorical_results:
            nonhistorical_results.append(variable_name)
        result_file_configuration["nodal_nonhistorical_results"].SetStringArray(nonhistorical_results)

        KM.Logger.PrintWarning(
            "GiDDamOutputProcess",
            "'NODAL_YOUNG_MODULUS' is not a historical variable. It is written as a "
            "nonhistorical nodal result prepared from element material properties."
        )

    def _translate_parameters(self, model_part, start_time, legacy_param):
        result_file_configuration = legacy_param["result_file_configuration"]
        output_control_type = result_file_configuration["output_control_type"].GetString()

        if output_control_type == "step":
            # The legacy implementation treated "step" as time in seconds.
            KM.Logger.PrintWarning(
                "GiDDamOutputProcess",
                "Legacy 'step' output is interpreted as time in seconds. "
                "Use 'time_s' for new configurations."
            )
            time_unit_factor = 1.0
        elif output_control_type in self._time_unit_factors:
            time_unit_factor = self._time_unit_factors[output_control_type]
        else:
            raise Exception(
                "GiDDamOutputProcess Error: Unknown value \"{}\" read for parameter "
                "\"output_control_type\"".format(output_control_type)
            )

        self._dam_output_factor = time_unit_factor
        self._dam_output_interval = result_file_configuration["output_interval"].GetDouble()
        self._dam_start_output_results = result_file_configuration["start_output_results"].GetDouble()
        # The legacy process mixed simulation and output units here. Convert the
        # simulation start time to output units; configurations whose units
        # coincide retain their previous behavior.
        self._dam_start_time_output = (
            start_time * self._get_simulation_time_factor(model_part) / time_unit_factor
        )

        result_file_configuration["output_control_type"].SetString("time")
        result_file_configuration["output_interval"].SetDouble(
            self._dam_output_interval * time_unit_factor)
        result_file_configuration.RemoveValue("start_output_results")
        result_file_configuration.AddEmptyValue("file_label").SetString("time")
        result_file_configuration.AddEmptyValue("time_label_format").SetString("{:.12f}")

        return legacy_param

    @staticmethod
    def _get_simulation_time_factor(model_part):
        try:
            from KratosMultiphysics.DamApplication import TIME_UNIT_CONVERTER
            return model_part.ProcessInfo[TIME_UNIT_CONVERTER]
        except RuntimeError:
            return 1.0

    def ExecuteInitialize(self):
        super().ExecuteInitialize()

        if self._dam_start_time_output >= self._dam_start_output_results:
            self._dam_next_output = (
                self._dam_start_time_output + self._dam_output_interval
            )
        else:
            self._dam_next_output = self._dam_start_output_results

        if self.model_part.ProcessInfo[KM.IS_RESTARTED]:
            self._dam_next_output = self._current_output_time()

    def ExecuteBeforeSolutionLoop(self):
        if self.multifile_flag == KM.MultiFileFlag.SingleFile:
            super().ExecuteBeforeSolutionLoop()
            if self._has_initial_result_frame():
                self._write_initial_result_values(self._initial_output_label())
        else:
            if self._has_initial_result_frame():
                self._write_initial_multifile_frame(self._initial_output_label())
            if self.point_output_process is not None:
                self.point_output_process.ExecuteBeforeSolutionLoop()

    def IsOutputStep(self):
        return self._current_output_time() >= self._dam_next_output

    def PrintOutput(self):
        super().PrintOutput()

        if self._dam_output_interval > 0.0:
            current_output_time = self._current_output_time()
            while self._dam_next_output <= current_output_time:
                self._dam_next_output += self._dam_output_interval

    def _current_output_time(self):
        return self.model_part.ProcessInfo[KM.TIME] / self._dam_output_factor

    def _has_initial_result_frame(self):
        return (
            self._dam_start_output_results == 0.0 or
            self._dam_start_time_output != 0.0
        )

    def _initial_output_label(self):
        return max(self._dam_start_output_results, self._dam_start_time_output)

    # The standard writer exposes these result writers only as private methods.
    # Reusing them avoids duplicating the writer; the adapter tests cover them.
    def _write_initial_result_values(self, label):
        self._GiDOutputProcess__write_nodal_results(label)
        self._GiDOutputProcess__write_gp_results(label)
        self._GiDOutputProcess__write_nonhistorical_nodal_results(label)
        self._GiDOutputProcess__write_nodal_flags(label)

    def _write_initial_multifile_frame(self, label):
        self._GiDOutputProcess__write_mesh(label)
        self._GiDOutputProcess__initialize_results(label)
        self._GiDOutputProcess__write_nodal_results(label)
        # The legacy MultipleFiles initial frame did not write Gauss-point results.
        self._GiDOutputProcess__write_nonhistorical_nodal_results(label)
        self._GiDOutputProcess__write_nodal_flags(label)
        self._GiDOutputProcess__finalize_results()
        self._GiDOutputProcess__write_step_to_list(label)

    def _get_legacy_pretty_time(self, time):
        return time / self._dam_output_factor
