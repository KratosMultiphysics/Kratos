import KratosMultiphysics as KM
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving
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

    The angular velocity may be applied either directly or after a linear ramp-up
    defined by an acceleration time.
    """

    def __init__(self, model, settings):
        KM.Process.__init__(self)

        default_settings = KM.Parameters("""{
            "rotating_frame_model_part_name": "",
            "rotating_object_model_part_name": "",
            "center_of_rotation": [0.0, 0.0, 0.0],
            "axis_of_rotation": [1.0, 0.0, 0.0],
            "target_angular_velocity_radians": 0.0,
            "acceleration_time": 0.0,
            "implementation": "cpp",
            "echo_level": 0,
            "fix_mesh_displacement": false,
            "fix_velocity": false
        }""")

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

        self.target_angular_velocity_radians = settings["target_angular_velocity_radians"].GetDouble()

        self.acceleration_time = settings["acceleration_time"].GetDouble()
        if self.acceleration_time < 0.0:
            raise Exception("The 'acceleration_time' parameter must be non-negative.")

        self.implementation = settings["implementation"].GetString().strip().lower()
        if self.implementation in ("c++", "cxx"):
            self.implementation = "cpp"
        if self.implementation not in ("cpp", "python"):
            raise Exception(
                f"Unsupported 'implementation': '{self.implementation}'. Use 'cpp' or 'python'."
            )

        self.echo_level = settings["echo_level"].GetInt()
        if self.echo_level < 0:
            raise Exception("The 'echo_level' parameter must be >= 0.")

        self.fix_mesh_displacement = settings["fix_mesh_displacement"].GetBool()
        self.fix_velocity = settings["fix_velocity"].GetBool()
        self._python_cache_initialized = False

    def ExecuteInitializeSolutionStep(self):
        current_time = self.rotating_frame_model_part.ProcessInfo[KM.TIME]

        theta, omega = self._GetThetaAndOmega(current_time)

        if self.implementation == "cpp":
            t0 = time.perf_counter()
            KratosMeshMoving.RotatingFrameUtility.ApplyRotationAndMeshDisplacement(
                self.rotating_frame_model_part,
                self.axis_of_rotation,
                theta,
                self.center_of_rotation
            )
            rot_time = time.perf_counter() - t0

            t0 = time.perf_counter()
            KratosMeshMoving.RotatingFrameUtility.AssignRotationalVelocity(
                self.rotating_object_model_part,
                self.axis_of_rotation,
                omega,
                self.center_of_rotation
            )
            vel_time = time.perf_counter() - t0
        else:
            if not self._python_cache_initialized:
                self._InitializePythonBackendCaches()

            t0 = time.perf_counter()
            self._ApplyRotationAndMeshDisplacementPython(theta)
            rot_time = time.perf_counter() - t0

            t0 = time.perf_counter()
            self._AssignRotationalVelocityPython(omega)
            vel_time = time.perf_counter() - t0

        if self.echo_level > 0:
            KM.Logger.PrintInfo(
                "RotatingFrameProcess",
                f"[implementation={self.implementation}] "
                f"rotation={rot_time*1.0e3:.3f} ms, "
                f"velocity={vel_time*1.0e3:.3f} ms"
            )

        if self.fix_mesh_displacement:
            self._FixVectorVariable(self.rotating_frame_model_part, KM.MESH_DISPLACEMENT)

        if self.fix_velocity:
            self._FixVectorVariable(self.rotating_object_model_part, KM.VELOCITY)

    def _GetThetaAndOmega(self, time):
        if np.isclose(self.acceleration_time, 0.0):
            omega = self.target_angular_velocity_radians
            theta = omega * time
            return theta, omega

        alpha = self.target_angular_velocity_radians / self.acceleration_time

        if time <= self.acceleration_time:
            omega = alpha * time
            theta = 0.5 * alpha * time**2
        else:
            omega = self.target_angular_velocity_radians
            theta = 0.5 * alpha * self.acceleration_time**2 + omega * (time - self.acceleration_time)

        return theta, omega

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
    def _FixVectorVariable(model_part, variable):
        components_map = {
            KM.VELOCITY: (KM.VELOCITY_X, KM.VELOCITY_Y, KM.VELOCITY_Z),
            KM.MESH_DISPLACEMENT: (KM.MESH_DISPLACEMENT_X, KM.MESH_DISPLACEMENT_Y, KM.MESH_DISPLACEMENT_Z),
        }

        if variable not in components_map:
            raise Exception(f"Fixing variable '{variable.Name()}' is not supported by this process.")

        vu = KM.VariableUtils()
        for comp in components_map[variable]:
            vu.ApplyFixity(comp, True, model_part.Nodes)

    def _ApplyRotationAndMeshDisplacementPython(self, theta):
        sin_half_theta = np.sin(theta / 2.0)
        a = np.cos(theta / 2.0)
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

        current_positions = self._rf_current_pos_ta.data
        np.subtract(self._rf_initial_positions, self.center_of_rotation, out=current_positions)
        current_positions[:] = current_positions @ rot_matrix
        current_positions += self.center_of_rotation

        np.subtract(current_positions, self._rf_initial_positions, out=self._rf_mesh_disp_ta.data)

        self._rf_mesh_disp_ta.StoreData()
        self._rf_current_pos_ta.StoreData()

    def _AssignRotationalVelocityPython(self, omega):
        angular_velocity_vector = omega * self.axis_of_rotation

        self._ro_current_pos_ta.CollectData()
        np.subtract(
            self._ro_current_pos_ta.data,
            self.center_of_rotation,
            out=self._ro_relative_positions
        )

        # v = omega x r, matching C++ MathUtils<double>::CrossProduct(omega, r)
        v = self._ro_velocity_ta.data
        r = self._ro_relative_positions
        v[:, 0] = angular_velocity_vector[1] * r[:, 2] - angular_velocity_vector[2] * r[:, 1]
        v[:, 1] = angular_velocity_vector[2] * r[:, 0] - angular_velocity_vector[0] * r[:, 2]
        v[:, 2] = angular_velocity_vector[0] * r[:, 1] - angular_velocity_vector[1] * r[:, 0]

        self._ro_velocity_ta.StoreData()
