from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.HDF5Application as KratosHDF5

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return HDF5OutputProcess(Model, settings["Parameters"])

# All the processes python processes should be derived from "Process"

class HDF5OutputProcess(KratosMultiphysics.Process):
    """This class is used in order to compute some pre and post process on the SPRISM solid shell elements

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing the settings.
    """

    def __init__(self, Model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """

        # Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "help"                              : "This process creates HDF5 output files",
            "model_part_name"                   : "",
            "create_xdmf_file"                  : true,
            "save_restart_files_in_folder"      : true,  # this should be named differently!
            "file_settings" : {
                "file_access_mode"              : "truncate"
            },
            "nodal_solution_step_data_settings" : {
                "list_of_variables": [ ]
            },
            "element_data_value_settings"       : {
                "list_of_variables": [ ]
            },
            "output_time_settings"              : {
                "output_step_frequency": 1
            }
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(default_parameters)

        # We define the model parts
        self.model_part = Model[self.settings["model_part_name"].GetString()]

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass
        # Create the actual output-processes
        is_mpi_execution = (self.model_part.GetCommunicator().TotalProcesses() > 1)

        self.hfd5_writer_process = ...

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteBeforeOutputStep(self):
        """ This method is executed right before the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def IsOutputStep(self):
        return True
        # if self.output_control_is_time:
        #     time = self.__get_pretty_time(self.model_part.ProcessInfo[TIME])
        #     return (time >= self.__get_pretty_time(self.next_output))
        # else:
        #     return ( self.step_count >= self.next_output )

    def PrintOutput(self):
        # here we have to call the ExecuteFinalizeSolutionStep of the processes bcs thats what they do
        pass

    def ExecuteAfterOutputStep(self):
        """ This method is executed right after the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteFinalize(self):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

        # Create xdmf-file
        if self.settings["create_xdmf_file"].GetBool():
            from create_xdmf_file import Execute # todo this method does not exist yet!
            Execute(file_name)
