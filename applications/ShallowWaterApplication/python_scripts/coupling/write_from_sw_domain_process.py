import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
from KratosMultiphysics.HDF5Application import single_mesh_temporal_output_process

def Factory(settings: KM.Parameters, model: KM.Model):
    if not isinstance(model, KM.Model):
        raise Exception("expected input shall be a Model object, encapsulating a json string")
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return WriteFromSwDomainProcess(model, settings["Parameters"])

class WriteFromSwDomainProcess(KM.OutputProcess):
    """WriteFromSwDomainProcess

    This process writes values from the shallow water domain to an HDF5 file.
    """

    @staticmethod
    def GetDefaultParameters():
        return KM.Parameters("""{
            "sw_model_part_name"            : "",
            "store_historical_database"     : true,
            "store_non_historical_database" : true,
            "interval"                      : [0.0,"End"],
            "output_control_type"           : "step",
            "output_interval"               : 1.0,
            "file_settings"                 : {
                "save_output_files_in_folder"   : true,
                "output_path"                   : "hdf5_output",
                "max_files_to_keep"             : "unlimited",
                "echo_level"                    : 0
            }
        }""")

    def __init__(self, model, settings):
        """The constructor of the WriteFromSwDomainProcess"""

        KM.OutputProcess.__init__(self)
        self.model = model
        print("hehe")
        print(self.model)
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.interval = KM.IntervalUtility(self.settings)

        self.hdf5_process = single_mesh_temporal_output_process.Factory(self._CreateHDF5Parameters(), model)


    def Check(self):
        '''Check the processes.'''
        self.hdf5_process.Check()


    def ExecuteBeforeSolutionLoop(self):
        '''Write the model part in HDF5 format.'''
        self.hdf5_process.ExecuteBeforeSolutionLoop()

        # this is not working for some reason
        # self.__controller = KM.OutputController(self.model, self.settings)


    def IsOutputStep(self):
        # return self.__controller.Evaluate() # this is not working for some reason
        return True


    def PrintOutput(self):
        '''Perform the depth integration over the interface model part.'''
        self.hdf5_process.ExecuteFinalizeSolutionStep()
        # this is not working for some reason
        # self.__controller.Update()


    def _CreateHDF5Parameters(self):
        hdf5_settings = KM.Parameters()
        hdf5_settings.AddValue("model_part_name", self.settings["sw_model_part_name"])
        hdf5_settings.AddValue("file_settings", self.settings["file_settings"])

        hdf5_settings.AddEmptyValue("output_time_settings")

        if self.settings["output_control_type"].GetString() == "step":
            output_interval = self.settings["output_interval"]
            hdf5_settings["output_time_settings"].AddValue("step_frequency", output_interval)

        elif self.settings["output_control_type"].GetString() == "time":
            output_interval = self.settings["output_interval"]
            hdf5_settings["output_time_settings"].AddValue("time_frequency", output_interval)


        historical_data_settings = KM.Parameters("""
        {"list_of_variables"       : ["MOMENTUM",
                                      "VELOCITY",
                                      "HEIGHT",
                                      "VERTICAL_VELOCITY",
                                      "TOPOGRAPHY"]}""")

        non_historical_data_settings = KM.Parameters("""{"list_of_variables": ["NODAL_AREA"]}""")

        if self.settings["store_historical_database"].GetBool():
            hdf5_settings.AddValue("nodal_solution_step_data_settings", historical_data_settings)
        # if self.settings["store_non_historical_database"].GetBool():
        #     hdf5_settings.AddValue("nodal_data_value_settings", non_historical_data_settings)

        hdf5_process_settings = KM.Parameters()
        hdf5_process_settings.AddValue("Parameters", hdf5_settings)

        return hdf5_process_settings
