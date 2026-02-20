# Importing the Kratos Library
import KratosMultiphysics

# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

try:
    from KratosMultiphysics.FluidDynamicsApplication import FlowForcesAndMomentsUtilities
except ImportError:
    FlowForcesAndMomentsUtilities = KratosMultiphysics.FlowForcesAndMomentsUtilities


def Factory(settings, model):
    if type(settings) != KratosMultiphysics.Parameters:
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ComputeAerodynamicCoefficientsProcess(model, settings["Parameters"])


class ComputeAerodynamicCoefficientsProcess(KratosMultiphysics.Process):
    """
    Compute aerodynamic force and moment coefficients from either:
      - body-fitted forces (REACTION-based), or
      - embedded forces (DRAG_FORCE-based),
    selected via the parameter "is_embedded".

    Coefficients:
        C_F = F / (q * Sref)
        C_M = M / (q * Sref * cref)
    """

    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters(r"""
        {
            "model_part_name"                 : "",
            "interval"                        : [0.0, 1e30],
            "print_coefficients_to_screen"    : false,
            "print_format"                    : ".8f",
            "write_output_file"               : true,
            "output_file_settings"            : {},
            "is_embedded"                     : false,
            "reference_point"                 : [0.0, 0.0, 0.0],
            "reference_area"                  : 1.0,
            "reference_chord"                 : 1.0,
            "dynamic_pressure"                : 1.0
        }
        """)

        self.params = params

        # Detect "End" as a tag and replace it by a large number 
        if self.params.Has("interval"):
            if self.params["interval"][1].IsString():
                if self.params["interval"][1].GetString() == "End":
                    self.params["interval"][1].SetDouble(1e30)
                else:
                    raise Exception(
                        'The second value of interval can be "End" or a number, interval currently: '
                        + self.params["interval"].PrettyPrintJsonString()
                    )

        self.params.ValidateAndAssignDefaults(default_settings)

        self.format = self.params["print_format"].GetString()

        # Get ModelPart
        model_part_name = self.params["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception('No "model_part_name" was specified!')
        self.model_part = model[model_part_name]

        self.flow_forces__and_moments_utilities = FlowForcesAndMomentsUtilities()

    def ExecuteInitialize(self):
        self.interval = KratosMultiphysics.Vector(2)
        self.interval[0] = self.params["interval"][0].GetDouble()
        self.interval[1] = self.params["interval"][1].GetDouble()

        self.print_to_screen = self.params["print_coefficients_to_screen"].GetBool()
        self.write_output_file = self.params["write_output_file"].GetBool()

        self.is_embedded = self.params["is_embedded"].GetBool()

        # Moment reference point
        self.reference_point = KratosMultiphysics.Vector(3)
        self.reference_point[0] = self.params["reference_point"][0].GetDouble()
        self.reference_point[1] = self.params["reference_point"][1].GetDouble()
        self.reference_point[2] = self.params["reference_point"][2].GetDouble()

        # Coefficient normalization data
        self.reference_area = self.params["reference_area"].GetDouble()
        self.reference_chord = self.params["reference_chord"].GetDouble()
        self.dynamic_pressure = self.params["dynamic_pressure"].GetDouble()

        if self.reference_area <= 0.0:
            raise Exception('"reference_area" must be > 0.0')
        if self.reference_chord <= 0.0:
            raise Exception('"reference_chord" must be > 0.0')
        if self.dynamic_pressure <= 0.0:
            raise Exception('"dynamic_pressure" must be > 0.0')

        # Output file
        if self.write_output_file and self.model_part.GetCommunicator().MyPID() == 0:
            # Default filename
            output_file_name = self.params["model_part_name"].GetString() + "_aero_coeffs.dat"

            file_handler_params = KratosMultiphysics.Parameters(self.params["output_file_settings"])

            if file_handler_params.Has("file_name"):
                warn_msg  = 'Unexpected user-specified entry found in "output_file_settings": {"file_name": '
                warn_msg += '"' + file_handler_params["file_name"].GetString() + '"}\n'
                warn_msg += 'Using this specified file name instead of the default "' + output_file_name + '"'
                KratosMultiphysics.Logger.PrintWarning("ComputeAerodynamicCoefficientsProcess", warn_msg)
            else:
                file_handler_params.AddEmptyValue("file_name")
                file_handler_params["file_name"].SetString(output_file_name)

            file_header = self._GetFileHeader()
            self.output_file = TimeBasedAsciiFileWriterUtility(
                self.model_part, file_handler_params, file_header
            ).file
            self.output_file.flush()

    def ExecuteFinalizeSolutionStep(self):
        current_step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if not (((current_time >= self.interval[0]) and (current_time < self.interval[1])) and
                (current_step + 1 >= self.model_part.GetBufferSize())):
            return

        # compute dimensional forces & moments 
        flow_force, flow_moment = self._GetFlowForcesAndMoments()

        denom_f = self.dynamic_pressure * self.reference_area
        denom_m = denom_f * self.reference_chord

        c_force = [
            flow_force[0] / denom_f,
            flow_force[1] / denom_f,
            flow_force[2] / denom_f
        ]
        c_moment = [
            flow_moment[0] / denom_m,
            flow_moment[1] / denom_m,
            flow_moment[2] / denom_m
        ]

        if self.model_part.GetCommunicator().MyPID() != 0:
            return

        if self.print_to_screen:
            msg = (
                str(current_time)
                + " Cx: " + format(c_force[0], self.format)
                + " Cy: " + format(c_force[1], self.format)
                + " Cz: " + format(c_force[2], self.format)
                + " Cmx: " + format(c_moment[0], self.format)
                + " Cmy: " + format(c_moment[1], self.format)
                + " Cmz: " + format(c_moment[2], self.format)
            )
            self._PrintToScreen(msg)

        if self.write_output_file:
            self.output_file.write(
                str(current_time) + " "
                + format(c_force[0], self.format) + " "
                + format(c_force[1], self.format) + " "
                + format(c_force[2], self.format) + " "
                + format(c_moment[0], self.format) + " "
                + format(c_moment[1], self.format) + " "
                + format(c_moment[2], self.format) + "\n"
            )
            self.output_file.flush()

    def ExecuteFinalize(self):
        if self.write_output_file and self.model_part.GetCommunicator().MyPID() == 0:
            self.output_file.close()

    # ----------------- internals -----------------

    def _GetFileHeader(self):
        return "# HEADER\n# Time Cx Cy Cz Cmx Cmy Cmz\n"

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo("ComputeAerodynamicCoefficientsProcess", result_msg)

    def _GetFlowForcesAndMoments(self):
        ref = self.reference_point

        if self.is_embedded:
            return self.flow_forces__and_moments_utilities.CalculateEmbeddedFlowForcesAndMoments(
                self.model_part, ref
            )
        else:
            return self.flow_forces__and_moments_utilities.CalculateBodyFittedFlowForcesAndMoments(
                self.model_part, ref
            )