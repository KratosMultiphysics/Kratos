# Importing the Kratos Library
import KratosMultiphysics

# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ComputeFlowForcesAndMomentsProcess(model, settings["Parameters"])


class ComputeFlowForcesAndMomentsProcess(KratosMultiphysics.Process):
    """
    Auxiliary base class to output total flow forces and moments
    over obstacles in fluid dynamics problems.
    A derived class needs to be implemented to be able to use
    this functionality, as calling the base class alone is not enough.
    """
    def __init__(self, model, params ):
        """
        Auxiliary class to output total flow forces over obstacles
        in fluid dynamics problems for a body fitted model part.
        """
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"           : "",
                "interval"                  : [0.0, 1e30],
                "print_flow_forces_and_moments_to_screen"      : false,
                "print_format"              : ".8f",
                "write_flow_forces_and_moments_output_file"    : true,
                "output_file_settings": {},
                "reference_point"           : [0.0, 0.0, 0.0]
            }
            """)

        self.params = params
        # Detect "End" as a tag and replace it by a large number
        if(self.params.Has("interval")):
            if(self.params["interval"][1].IsString()):
                if(self.params["interval"][1].GetString() == "End"):
                    self.params["interval"][1].SetDouble(1e30)
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:"+self.params["interval"].PrettyPrintJsonString())

        self.params.ValidateAndAssignDefaults(default_settings)

        self.format = self.params["print_format"].GetString()

        # getting the ModelPart from the Model
        model_part_name = self.params["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception('No "model_part_name" was specified!')
        else:
            self.model_part = model[self.params["model_part_name"].GetString()]

    def ExecuteInitialize(self):

        self.interval = KratosMultiphysics.Vector(2)
        self.interval[0] = self.params["interval"][0].GetDouble()
        self.interval[1] = self.params["interval"][1].GetDouble()
        self.print_flow_forces_and_moments_to_screen = self.params["print_flow_forces_and_moments_to_screen"].GetBool()
        self.write_flow_forces_and_moments_output_file = self.params["write_flow_forces_and_moments_output_file"].GetBool()
        self.reference_point = KratosMultiphysics.Vector(3)
        self.reference_point[0] = self.params["reference_point"][0].GetDouble()
        self.reference_point[1] = self.params["reference_point"][1].GetDouble()
        self.reference_point[2] = self.params["reference_point"][2].GetDouble()

        if (self.write_flow_forces_and_moments_output_file):
            if (self.model_part.GetCommunicator().MyPID() == 0):

                output_file_name = self.params["model_part_name"].GetString() + "_drag.dat"

                file_handler_params = KratosMultiphysics.Parameters(self.params["output_file_settings"])

                if file_handler_params.Has("file_name"):
                    warn_msg  = 'Unexpected user-specified entry found in "output_file_settings": {"file_name": '
                    warn_msg += '"' + file_handler_params["file_name"].GetString() + '"}\n'
                    warn_msg += 'Using this specified file name instead of the default "' + output_file_name + '"'
                    KratosMultiphysics.Logger.PrintWarning("ComputeDragProcess", warn_msg)
                else:
                    file_handler_params.AddEmptyValue("file_name")
                    file_handler_params["file_name"].SetString(output_file_name)
                file_header = self._GetFileHeader()
                self.output_file = TimeBasedAsciiFileWriterUtility(self.model_part,
                    file_handler_params, file_header).file
                self.output_file.flush()

    def ExecuteFinalizeSolutionStep(self):
        current_step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if(((current_time >= self.interval[0]) and  (current_time < self.interval[1])) and (current_step + 1 >= self.model_part.GetBufferSize())):
            # Compute the drag force
            flow_force, flow_moment = self._GetCorrespondingFlowForcesAndMoments()

            # Write the flow force values
            if (self.model_part.GetCommunicator().MyPID() == 0):
                if (self.print_flow_forces_and_moments_to_screen):
                    result_msg = str(current_time) + " x-force: " + format(flow_force[0],self.format) + " y-force: " + format(flow_force[1],self.format) + " z-force: " + format(flow_force[2],self.format) + \
                    " Mx: " + format(flow_moment[0],self.format) + " My: " + format(flow_moment[1],self.format) + " Mz: " + format(flow_moment[2],self.format)
                    self._PrintToScreen(result_msg)

                # not formatting time in order to not lead to problems with time recognition
                # in the file writer when restarting
                if (self.write_flow_forces_and_moments_output_file):
                    self.output_file.write(str(current_time)+" "+format(flow_force[0],self.format)+" "+format(flow_force[1],self.format)+" "+format(flow_force[2],self.format)+" "+ \
                    format(flow_moment[0],self.format)+" "+format(flow_moment[1],self.format)+" "+format(flow_moment[2],self.format)+"\n")
                    self.output_file.flush()

    def ExecuteFinalize(self):
        if (self.write_flow_forces_and_moments_output_file):
            if (self.model_part.GetCommunicator().MyPID() == 0):
                self.output_file.close()

    def _GetFileHeader(self):
        err_msg  = 'ComputeFlowForcesAndMomentsProcess: _GetFileHeader called in base class\n'
        err_msg += 'this needs to be implemented and called from derived class'
        raise Exception(err_msg)

    def _PrintToScreen(self,result_msg):
        err_msg  = 'ComputeFlowForcesAndMomentsProcess: _PrinToScreen called in base class\n'
        err_msg += 'this needs to be implemented and called from derived class'
        raise Exception(err_msg)

    def _GetCorrespondingFlowForcesAndMoments(self):
        err_msg  = 'ComputeFlowForcesAndMomentsProcess: _GetCorrespondingFlowForcesAndMoments called in base class\n'
        err_msg += 'this needs to be implemented and called from derived class'
        raise Exception(err_msg)