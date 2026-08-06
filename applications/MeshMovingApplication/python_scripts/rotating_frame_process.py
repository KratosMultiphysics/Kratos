import KratosMultiphysics as KM
import numpy as np
import time


def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return RotatingFrameProcess(model, settings["Parameters"])


class RotatingFrameProcess(KM.Process):
    """Process to apply rigid-body rotational kinematics to model parts.

    This process:
      - rotates a model part and assigns the corresponding MESH_DISPLACEMENT
      - assigns the corresponding rigid-body rotational VELOCITY to another model part

    The angular velocity is provided as the standard scalar input
    'angular_velocity_radians = f(x,y,z,t)' (or a constant number).
    The rotation is advanced by integrating the angular-velocity increment
    at each solution step.
    """

    def __init__(self, model, settings):
        KM.Process.__init__(self)

        default_settings = KM.Parameters("""{
            "rotating_frame_model_part_name": "",
            "rotating_object_model_part_name": "",
            "interval": [0.0, 1e30],
            "center_of_rotation": [0.0, 0.0, 0.0],
            "axis_of_rotation": [1.0, 0.0, 0.0],
            "angular_velocity_radians": 0.0,
            "echo_level": 0,
            "fix_mesh_displacement": false,
            "fix_velocity": false
        }""")

        # Allow a string function expression for 'angular_velocity_radians'.
        if settings.Has("angular_velocity_radians") and settings["angular_velocity_radians"].IsString():
            default_settings["angular_velocity_radians"].SetString("0.0")

        # Assign this here since it may normalize the "interval" before validation
        self.interval = KM.IntervalUtility(settings)

        settings.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.settings = settings

        if not settings["rotating_frame_model_part_name"].GetString():
            raise Exception("'rotating_frame_model_part_name' not provided. Please specify the model part to be rotated.")
        self.rotating_frame_model_part_name = settings["rotating_frame_model_part_name"].GetString()

        if not settings["rotating_object_model_part_name"].GetString():
            raise Exception("'rotating_object_model_part_name' not provided. Please specify the model part to which rotational velocity will be assigned.")
        self.rotating_object_model_part_name = settings["rotating_object_model_part_name"].GetString()

        self.rotating_frame_model_part = self.model.GetModelPart(self.rotating_frame_model_part_name)
        self.rotating_object_model_part = self.model.GetModelPart(self.rotating_object_model_part_name)
        self._interval_begin = settings["interval"][0].GetDouble()
        self._interval_end = settings["interval"][1].GetDouble()

        self.center_of_rotation = np.array(settings["center_of_rotation"].GetVector(), dtype=float)

        self.axis_of_rotation = np.array(settings["axis_of_rotation"].GetVector(), dtype=float)
        if self.axis_of_rotation.size == 0:
            raise Exception("The 'axis_of_rotation' vector is empty.")

        axis_norm = np.linalg.norm(self.axis_of_rotation)
        if np.isclose(axis_norm, 0.0):
            raise Exception("The 'axis_of_rotation' vector must have non-zero norm.")

        if not np.isclose(axis_norm, 1.0, rtol=1e-6):
            KM.Logger.PrintWarning(
                "RotatingFrameProcess",
                "The 'axis_of_rotation' vector is not a unit vector. It will be normalized."
            )
            self.axis_of_rotation /= axis_norm

        self.angular_velocity_function = self._CreateScalarFunction(
            settings["angular_velocity_radians"],
            "angular_velocity_radians"
        )

        self.echo_level = settings["echo_level"].GetInt()
        if self.echo_level < 0:
            raise Exception("The 'echo_level' parameter must be >= 0.")

        self.fix_mesh_displacement = settings["fix_mesh_displacement"].GetBool()
        self.fix_velocity = settings["fix_velocity"].GetBool()
        self._fixity_applied = False
        self._python_cache_initialized = False

    def ExecuteInitializeSolutionStep(self):
        current_time = self.rotating_frame_model_part.ProcessInfo[KM.TIME]
        delta_theta, omega, is_active = self._GetDeltaThetaAndOmegaFromFunction(current_time)

        if np.isclose(delta_theta, 0.0) and not is_active:
            # If we leave the active interval, release any fixity that this process applied.
            if self._fixity_applied:
                if self.fix_mesh_displacement:
                    self._FixVectorVariable(self.rotating_frame_model_part, KM.MESH_DISPLACEMENT, False)
                if self.fix_velocity:
                    self._FixVectorVariable(self.rotating_object_model_part, KM.VELOCITY, False)
                self._fixity_applied = False
            return

        if not self._python_cache_initialized:
            self._InitializePythonBackendCaches()

        t0 = time.perf_counter()
        if not np.isclose(delta_theta, 0.0):
            self._ApplyRotationAndMeshDisplacementPython(delta_theta)
        rot_time = time.perf_counter() - t0

        t0 = time.perf_counter()
        if is_active:
            self._AssignRotationalVelocityPython(omega)
        vel_time = time.perf_counter() - t0

        if self.echo_level > 0:
            KM.Logger.PrintInfo(
                "RotatingFrameProcess", f" rotation={rot_time*1.0e3:.3f} ms, velocity={vel_time*1.0e3:.3f} ms"
            )

        if is_active:
            if self.fix_mesh_displacement:
                self._FixVectorVariable(self.rotating_frame_model_part, KM.MESH_DISPLACEMENT, True)

            if self.fix_velocity:
                self._FixVectorVariable(self.rotating_object_model_part, KM.VELOCITY, True)

            self._fixity_applied = self.fix_mesh_displacement or self.fix_velocity
        elif self._fixity_applied:
            if self.fix_mesh_displacement:
                self._FixVectorVariable(self.rotating_frame_model_part, KM.MESH_DISPLACEMENT, False)
            if self.fix_velocity:
                self._FixVectorVariable(self.rotating_object_model_part, KM.VELOCITY, False)
            self._fixity_applied = False

    @staticmethod
    def _CreateScalarFunction(value, value_name):
        if value.IsNumber():
            function_string = str(value.GetDouble())
        elif value.IsString():
            function_string = value.GetString()
            if function_string == "":
                raise Exception(f"The '{value_name}' function string is empty.")
        else:
            raise Exception(
                f"The '{value_name}' parameter must be provided as a number or as a function string."
            )
        return KM.GenericFunctionUtility(function_string)

    def _EvaluateAngularVelocityFunction(self, time):
        x = self.center_of_rotation[0]
        y = self.center_of_rotation[1]
        z = self.center_of_rotation[2]
        return self.angular_velocity_function.CallFunction(x, y, z, time, 0.0, 0.0, 0.0)

    def _GetDeltaThetaAndOmegaFromFunction(self, time):
        omega = self._EvaluateAngularVelocityFunction(time)
        delta_time = self.rotating_frame_model_part.ProcessInfo[KM.DELTA_TIME]
        delta_theta = 0.0
        if delta_time > 0.0:
            step_begin = time - delta_time
            overlap_begin = max(step_begin, self._interval_begin)
            overlap_end = min(time, self._interval_end)

            if overlap_end > overlap_begin:
                omega_begin = self._EvaluateAngularVelocityFunction(overlap_begin)
                omega_end = self._EvaluateAngularVelocityFunction(overlap_end)
                delta_theta = 0.5 * (omega_begin + omega_end) * (overlap_end - overlap_begin)

        return delta_theta, omega, self.interval.IsInInterval(time)

    def _InitializePythonBackendCaches(self):
        self._rf_initial_pos_ta = KM.TensorAdaptors.NodePositionTensorAdaptor(
            self.rotating_frame_model_part.Nodes, KM.Configuration.Initial
        )
        self._rf_current_pos_ta = KM.TensorAdaptors.NodePositionTensorAdaptor(
            self.rotating_frame_model_part.Nodes, KM.Configuration.Current
        )
        self._rf_mesh_disp_ta = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(
            self.rotating_frame_model_part.Nodes, KM.MESH_DISPLACEMENT
        )
        self._ro_current_pos_ta = KM.TensorAdaptors.NodePositionTensorAdaptor(
            self.rotating_object_model_part.Nodes, KM.Configuration.Current
        )
        self._ro_velocity_ta = KM.TensorAdaptors.HistoricalVariableTensorAdaptor(
            self.rotating_object_model_part.Nodes, KM.VELOCITY
        )

        # Initial frame coordinates are constant for a static mesh.
        self._rf_initial_pos_ta.CollectData()
        self._rf_initial_positions = self._rf_initial_pos_ta.data.copy()

        self._rf_rot_matrix = np.empty((3, 3), dtype=float)
        self._ro_relative_positions = np.empty(
            (self.rotating_object_model_part.NumberOfNodes(), 3),
            dtype=float
        )

        self._python_cache_initialized = True

    @staticmethod
    def _FixVectorVariable(model_part, variable, is_fixed):
        components_map = {
            KM.VELOCITY: (KM.VELOCITY_X, KM.VELOCITY_Y, KM.VELOCITY_Z),
            KM.MESH_DISPLACEMENT: (KM.MESH_DISPLACEMENT_X, KM.MESH_DISPLACEMENT_Y, KM.MESH_DISPLACEMENT_Z),
        }

        if variable not in components_map:
            raise Exception(f"Fixing variable '{variable.Name()}' is not supported by this process.")

        vu = KM.VariableUtils()
        for comp in components_map[variable]:
            vu.ApplyFixity(comp, is_fixed, model_part.Nodes)

    def _ApplyRotationAndMeshDisplacementPython(self, delta_theta):
        sin_half_theta = np.sin(delta_theta / 2.0)
        a = np.cos(delta_theta / 2.0)
        b = -self.axis_of_rotation[0] * sin_half_theta
        c = -self.axis_of_rotation[1] * sin_half_theta
        d = -self.axis_of_rotation[2] * sin_half_theta

        rot_matrix = self._rf_rot_matrix
        rot_matrix[0, 0] = a * a + b * b - c * c - d * d
        rot_matrix[0, 1] = 2.0 * (b * c - a * d)
        rot_matrix[0, 2] = 2.0 * (b * d + a * c)
        rot_matrix[1, 0] = 2.0 * (b * c + a * d)
        rot_matrix[1, 1] = a * a + c * c - b * b - d * d
        rot_matrix[1, 2] = 2.0 * (c * d - a * b)
        rot_matrix[2, 0] = 2.0 * (b * d - a * c)
        rot_matrix[2, 1] = 2.0 * (c * d + a * b)
        rot_matrix[2, 2] = a * a + d * d - b * b - c * c

        # Rotate current coordinates around the center with the incremental angle.
        self._rf_current_pos_ta.CollectData()
        current_positions = self._rf_current_pos_ta.data
        current_positions[:] = (current_positions - self.center_of_rotation) @ rot_matrix
        current_positions += self.center_of_rotation

        # Mesh displacement from reference configuration: u_mesh = x - x0
        self._rf_mesh_disp_ta.data[:] = current_positions - self._rf_initial_positions

        self._rf_mesh_disp_ta.StoreData()
        self._rf_current_pos_ta.StoreData()

    def _AssignRotationalVelocityPython(self, omega):
        angular_velocity_vector = omega * self.axis_of_rotation

        self._ro_current_pos_ta.CollectData()
        self._ro_relative_positions[:] = self._ro_current_pos_ta.data - self.center_of_rotation

        # v = omega x r, matching C++ MathUtils<double>::CrossProduct(omega, r)
        v = self._ro_velocity_ta.data
        r = self._ro_relative_positions
        v[:, 0] = angular_velocity_vector[1] * r[:, 2] - angular_velocity_vector[2] * r[:, 1]
        v[:, 1] = angular_velocity_vector[2] * r[:, 0] - angular_velocity_vector[0] * r[:, 2]
        v[:, 2] = angular_velocity_vector[0] * r[:, 1] - angular_velocity_vector[1] * r[:, 0]

        self._ro_velocity_ta.StoreData()
