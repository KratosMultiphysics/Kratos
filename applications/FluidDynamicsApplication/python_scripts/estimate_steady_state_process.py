# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# other imports
from time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return EstimateSteadyStateProcess(model, settings["Parameters"])


class EstimateSteadyStateProcess(KratosMultiphysics.Process):
    """
    Auxiliary base class to output total flow forces
    over obstacles in fluid dynamics problems.
    A derived class needs to be implemented to be able to use
    this functionality, as calling the base class alone is not enough.
    """
    def __init__(self, model, settings ):
        """
        Auxiliary class to output total flow forces over obstacles
        in fluid dynamics problems for a body fitted model part.
        """
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"           : "",
                "print_change_to_screen"     : false,
                "write_change_output_file"    : true,
                "print_format"              : ".8f",
                "output_file_settings": {}
            }
            """)

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.format = self.settings["print_format"].GetString()

        # getting the ModelPart from the Model
        model_part_name = self.settings["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception('No "model_part_name" was specified!')
        else:
            self.model_part = model[self.settings["model_part_name"].GetString()]

    def ExecuteInitialize(self):
        self.domain_size = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        # Compute the fluid domain NODAL_AREA values (required as weight for steady state estimation)
        KratosMultiphysics.CalculateNodalAreaProcess(self.model_part, 
                                                     self.domain_size).Execute()

        self.print_change_to_screen = self.settings["print_change_to_screen"].GetBool()
        self.write_change_output_file = self.settings["write_change_output_file"].GetBool()

        if (self.model_part.GetCommunicator().MyPID() == 0):
            if (self.write_change_output_file):

                output_file_name = self.settings["model_part_name"].GetString() + "_quantity_variation.dat"

                file_handler_settings = KratosMultiphysics.Parameters(self.settings["output_file_settings"])

                if file_handler_settings.Has("file_name"):
                    warn_msg  = 'Unexpected user-specified entry found in "output_file_settings": {"file_name": '
                    warn_msg += '"' + file_handler_settings["file_name"].GetString() + '"}\n'
                    warn_msg += 'Using this specififed file name instead of the default "' + output_file_name + '"'
                    KratosMultiphysics.Logger.PrintWarning("EstimateSteadyStateProcess", warn_msg)
                else:
                    file_handler_settings.AddEmptyValue("file_name")
                    file_handler_settings["file_name"].SetString(output_file_name)
                file_header = self._GetFileHeader()
                self.output_file = TimeBasedAsciiFileWriterUtility(self.model_part,
                    file_handler_settings, file_header).file

    def ExecuteFinalizeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        # Initializing steady state indicator utility
        SteadyStateIndicatorUtility = KratosCFD.SteadyStateIndicatorUtility(self.model_part)
        # Compute the quantity changes
        SteadyStateIndicatorUtility.EstimateQuantityChangesInTime()
        change_in_velocity = SteadyStateIndicatorUtility.GetVelocityChange()
        change_in_pressure = SteadyStateIndicatorUtility.GetPressureChange()

        # Write the quantity changes
        if (self.model_part.GetCommunicator().MyPID() == 0):
            if (self.print_change_to_screen):
                result_msg = str(current_time)
                result_msg +="\nChange in velocity in percentage: " + str(change_in_velocity) + "\n"
                result_msg +=  "Change in pressure in percentage: " + str(change_in_pressure) + "\n"
                self._PrintToScreen(result_msg)

            # not formatting time in order to not lead to problems with time recognition
            # in the file writer when restarting
            if (self.write_change_output_file):
                self.output_file.write(" change in velocity magnitude: " + format(change_in_velocity,self.format) + " change in pressure: " + format(change_in_pressure,self.format)+"\n")

    def ExecuteFinalize(self):
        if (self.model_part.GetCommunicator().MyPID() == 0):
            self.output_file.close()

    def _GetFileHeader(self):
        header  = '# Time Averaged Quantity Changes for Model Part ' + self.settings["model_part_name"].GetString() + '\n'
        header += '# Time Time-Averaged-Velocity-Magnitude Time-Averaged-Pressure\n'
        return header

    def _PrintToScreen(self,result_msg):
        KratosMultiphysics.Logger.PrintInfo("EstimateSteadyStateProcess","Time Averaged Quantity Change Results:")
        KratosMultiphysics.Logger.PrintInfo("EstimateSteadyStateProcess","\nCurrent time: " + result_msg)
