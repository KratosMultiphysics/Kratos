import KratosMultiphysics
import KratosMultiphysics.MeshMovingApplication as MeshMoving
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

import math
import sys



def Factory(parameters, model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise TypeError("expecting a Parameters object encapsulating a json string, but got {}".format(type(parameters)))
    return SetMeshMotionAndGetForcesProcess(model, parameters["Parameters"])



class SetMeshMotionAndGetForcesProcess(KratosMultiphysics.Process):
    def __init__(self, model, parameters):
        KratosMultiphysics.Process.__init__(self)

        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model_part = model[self.parameters["model_part_name"].GetString()]

        # Transformation variables
        self.reference_center = self.parameters["reference_point"].GetVector()
        self.reference_theta = self.parameters["z_rotation_angle"].GetDouble()
        self.updated_theta = self.reference_theta

        # Mesh motion parameters
        self.prescribed_motion = SetMeshMotionAndGetForcesProcess.Motion()
        self.prescribed_motion.pitch.amplitudes = [math.radians(value) for value in self.parameters["imposed_motion"]["pitch"]["amplitude"].GetVector()]
        self.prescribed_motion.pitch.circular_frequencies = [2.0 * math.pi * value for value in self.parameters["imposed_motion"]["pitch"]["frequency"].GetVector()]
        self.prescribed_motion.heave.amplitudes = [math.radians(value) for value in self.parameters["imposed_motion"]["heave"]["amplitude"].GetVector()]
        self.prescribed_motion.heave.circular_frequencies = [2.0 * math.pi * value for value in self.parameters["imposed_motion"]["heave"]["frequency"].GetVector()]
        self.prescribed_motion.rampup_time = parameters['rampup_time'].GetDouble()

        # Output parameters
        self.print_to_screen = self.parameters['print_to_screen'].GetBool()
        self.write_output_file = self.parameters['write_output_file'].GetBool()
        self.format = self.parameters["print_format"].GetString()
        self.rampup_time = self.parameters['rampup_time'].GetDouble()
        self.forced_flush_step_frequency = self.parameters["forced_flush_step_frequency"].GetInt()
        self.output_file = {}

        if self.parameters["output_file_settings"].Has("write_buffer_size"):
            if self.parameters["output_file_settings"]["write_buffer_size"].GetInt() < 0 and self.forced_flush_step_frequency > 0:
                self.parameters["output_file_settings"]["write_buffer_size"].SetInt(10)

        self._InitializeOutput(self.parameters)

        self.interval_utility = KratosMultiphysics.IntervalUtility(self.parameters)


    def ExecuteInitializeSolutionStep(self):
        # Evaluate imposed motion
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        pitch_increment = self.prescribed_motion.pitch.Evaluate(time)
        heave_increment = self.prescribed_motion.heave.Evaluate(time)

        rotation_axis = [0.0, 0.0, 1.0]
        rotation_angle = self.reference_theta + pitch_increment
        reference_point = self.reference_center
        translation_vector = [0.0, heave_increment, 0.0]

        # Execute mesh motion
        MeshMoving.MoveModelPart(
            self.model_part,
            rotation_axis,
            rotation_angle,
            reference_point,
            translation_vector)


    def ExecuteFinalizeSolutionStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]

        if self.interval_utility.IsInInterval(time):
            # Update reference point for moment reduction
            pitch_increment = self.prescribed_motion.pitch.Evaluate(time)
            heave_increment = self.prescribed_motion.heave.Evaluate(time)
            reference_point = KratosMultiphysics.Array3(self.reference_center) + KratosMultiphysics.Array3([0.0, heave_increment, 0.0])

            force, moment = KratosMultiphysics.ForceAndTorqueUtils.ComputeEquivalentForceAndTorque(
                self.model_part,
                reference_point,
                KratosMultiphysics.REACTION,
                KratosMultiphysics.MOMENT)

            force *= -1.0
            moment *= -1.0

            # Compute body-attached reactions
            rotation_axis = [0.0, 0.0, 1.0]
            rotation_angle = self.reference_theta + pitch_increment
            reference_point = [0.0, 0.0, 0.0]
            translation_vector = [0.0, 0.0, 0.0]

            transform = MeshMoving.AffineTransform(
                rotation_axis,
                rotation_angle,
                reference_point,
                translation_vector)

            body_force = transform.Apply(force)
            body_moment = transform.Apply(moment)

            # Write output
            if self.model_part.GetCommunicator().MyPID() == 0:
                output_values = list(force) + list(moment) + list(body_force) + list(body_moment)
                formatted_reactions = [format(value, self.format) for value in output_values]
                formatted_motion = [format(value, self.format) for value in [0.0, heave_increment, pitch_increment]]

                formatted_reactions.insert(0, str(time))
                formatted_motion.insert(0, str(time))

                if self.print_to_screen:
                    message = ["", ""]

                    message[0] += "Force evaluation for model part '{}'\n".format(self.model_part.Name)
                    message[0] += ", ".join([label+value for label, value in zip(self.reaction_labels, formatted_reactions)])

                    message[1] += "Prescribed motion for model part '{}'\n".format(self.model_part.Name)
                    message[1] += ", ".join([label + value for label, value in zip(self.motion_labels, formatted_motion)])

                    self._PrintToScreen(message)
                    sys.stdout.flush()

                if self.write_output_file:
                    self.output_file["force"].write(" ".join(formatted_reactions) + "\n")
                    self.output_file["motion"].write(" ".join(formatted_motion) + "\n")

                    if self.forced_flush_step_frequency > 0 and step % self.forced_flush_step_frequency == 0:
                        for key in self.output_file:
                            self.output_file[key].flush()


    @property
    def reaction_labels(self):
        return ['time: ',
                'fx: ', 'fy: ', 'fz: ', 'mx: ', 'my: ', 'mz: ',
                'fx\': ', 'fy\': ', 'fz\': ', 'mx\': ', 'my\': ', 'mz\': ']


    @property
    def motion_labels(self):
        return ['time: ',
                'dx: ', 'dy: ', 'dtheta: ']


    def _InitializeOutput(self, parameters):
        if self.write_output_file:
            output_file_name = parameters["model_part_name"].GetString()
            file_handler_parameters = KratosMultiphysics.Parameters(parameters["output_file_settings"])

            for case in ["motion", "force"]:
                case_file_handler_parameters = file_handler_parameters.Clone()

                if file_handler_parameters.Has("file_name"):
                    output_file_name = file_handler_parameters["file_name"].GetString()
                    warning_message = "Unexpected user-specified entry in 'output_file_settings': 'file_name':{}\n".format(file_handler_parameters["file_name"])
                    warning_message += "Using this specified file name instead of the default '{}'".format(output_file_name)
                    KratosMultiphysics.Logger.PrintWarning("SetMeshMotionAndWriteForcesProcess", warning_message)
                else:
                    case_file_handler_parameters.AddEmptyValue("file_name")

                case_file_handler_parameters["file_name"].SetString("{}_{}.dat".format(
                    output_file_name,
                    case))

                file_header = self._GetFileHeader(case)
                self.output_file[case] = TimeBasedAsciiFileWriterUtility(
                    self.model_part,
                    case_file_handler_parameters,
                    file_header).file


    def _GetFileHeader(self, case):
        header = {}

        header['force'] = '# Forces for model part ' + \
            self.model_part.Name + '\n'
        header['force'] += '# Time Fx Fy Fz Mx My Mz Fx\' Fy\' Fz\' Mx\' My\' Mz\'\n'

        header['motion'] = '# Motion for model part ' + \
            self.model_part.Name + '\n'
        header['motion'] += '# Time dx dy dtheta\n'

        return header[case]


    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo(
            'SetMeshMotionAndGetForcesProcess', 'Force - flow- and body-attached:')
        KratosMultiphysics.Logger.PrintInfo(
            'SetMeshMotionAndGetForcesProcess', 'Current time: ' + result_msg[0])
        KratosMultiphysics.Logger.PrintInfo(
            'SetMeshMotionAndGetForcesProcess', 'Motion:')
        KratosMultiphysics.Logger.PrintInfo(
            'SetMeshMotionAndGetForcesProcess', 'Current time: ' + result_msg[1])


    @staticmethod
    def GetDefaultParameters():
        return KratosMultiphysics.Parameters("""{
            "model_part_name" : "",
            "interval"              : [0.0, "End"],
            "rampup_time"           : 0.0,
            "reference_point"       : [0.0,0.0,0.0],
            "z_rotation_angle"      : 0.0,
            "imposed_motion":{
                "pitch": {"amplitude": [0.1], "frequency" : [1.2]},
                "heave": {"amplitude": [0.02, 0.03], "frequency" : [0.9, 1.1]}
            },
            "print_to_screen"       : false,
            "print_format"          : ".8f",
            "write_output_file"     : true,
            "output_file_settings"  : {},
            "forced_flush_step_frequency" : -1
        }""")


    class Motion:
        def __init__(self):
            self._rampup_time = 0.0
            self.pitch = SetMeshMotionAndGetForcesProcess.Motion.RampedSin([0.0], [0.0], self.rampup_time)
            self.heave = SetMeshMotionAndGetForcesProcess.Motion.RampedSin([0.0], [0.0], self.rampup_time)

        @property
        def rampup_time(self):
            return self._rampup_time

        @rampup_time.setter
        def rampup_time(self, value: float):
            self._rampup_time = value
            self.heave.rampup_time = self.rampup_time
            self.pitch.rampup_time = self.rampup_time

        class RampedSin:
            def __init__(self,
                        amplitudes: list,
                        circularFrequencies: list,
                        rampupTime: float):
                self.amplitudes = amplitudes
                self.circular_frequencies = circularFrequencies
                self.rampup_time = rampupTime

            def Evaluate(self, time: float):
                value = sum([amplitude * math.sin(circular_frequency * time) for amplitude, circular_frequency in zip(self.amplitudes, self.circular_frequencies)])
                if time <= self.rampup_time:
                    value *= 0.5 * (math.cos(math.pi * (1.0 - time/self.rampup_time)) + 1.0)
                return value