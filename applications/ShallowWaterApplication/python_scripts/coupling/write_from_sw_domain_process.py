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
            "output_model_part_name"        : "sw_output_model_part",
            "store_historical_database"     : true,
            "store_non_historical_database" : true,
            "interval"                      : [0.0,"End"],
            "output_control_type"           : "step",
            "output_interval"               : 1.0,
            "file_settings"                 : {
                "file_name": "hdf5_output/<model_part_name>-<time>.h5",
                "file_access_mode"              : "truncate",
                "max_files_to_keep"             : "unlimited",
                "echo_level"                    : 0
            }
        }""")

    def __init__(self, model, settings):
        """The constructor of the WriteFromSwDomainProcess"""

        KM.OutputProcess.__init__(self)
        self.model = model
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.variables = [KM.VELOCITY, KM.ACCELERATION, KM.MOMENTUM, SW.HEIGHT, SW.VERTICAL_VELOCITY, SW.TOPOGRAPHY]
        self.sw_model_part = model[self.settings["sw_model_part_name"].GetString()]
        self.output_model_part = model.CreateModelPart(self.settings["output_model_part_name"].GetString())

        self.interval = KM.IntervalUtility(self.settings)

        self.hdf5_process = single_mesh_temporal_output_process.Factory(self._CreateHDF5Parameters(), model)


    def Check(self):
        '''Check the processes.'''
        self.hdf5_process.Check()


    def ExecuteInitialize(self):
        '''Initialize the output model part.'''
        self._InitializeOutputModelPart()
        self._SetOutputProcessInfo()
        if not self.settings["store_historical_database"].GetBool():
            for var in self.variables:
                KM.VariableUtils().SetNonHistoricalVariableToZero(var, self.output_model_part.Nodes)


    def ExecuteBeforeSolutionLoop(self):
        '''Write the model part in HDF5 format.'''
        self._MapToOutputModelPart()
        self.hdf5_process.ExecuteBeforeSolutionLoop()
        # this is not working for some reason
        # self.__controller = KM.OutputController(self.model, self.settings)

    def ExecuteInitializeSolutionStep(self):
        '''Synchronize the ProcessInfo of the output and interface model part.'''
        self._SetOutputProcessInfo()

    def IsOutputStep(self):
        # return self.__controller.Evaluate() # this is not working for some reason
        return True


    def PrintOutput(self):
        '''Perform the depth integration over the interface model part.'''
        self._MapToOutputModelPart()
        self.hdf5_process.ExecuteFinalizeSolutionStep()
        # this is not working for some reason
        # self.__controller.Update()


    def _CreateHDF5Parameters(self):
        hdf5_settings = KM.Parameters()
        hdf5_settings.AddValue("model_part_name", self.settings["output_model_part_name"])
        hdf5_settings.AddValue("file_settings", self.settings["file_settings"])

        hdf5_settings.AddEmptyValue("output_time_settings")

        if self.settings["output_control_type"].GetString() == "step":
            output_interval = self.settings["output_interval"]
            hdf5_settings["output_time_settings"].AddValue("step_frequency", output_interval)

        elif self.settings["output_control_type"].GetString() == "time":
            output_interval = self.settings["output_interval"]
            hdf5_settings["output_time_settings"].AddValue("time_frequency", output_interval)


        historical_data_settings = KM.Parameters("""
        {"list_of_variables": [ "MOMENTUM",
                                "VELOCITY",
                                "ACCELERATION",
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

    def _InitializeOutputModelPart(self):
        if self.settings["store_historical_database"].GetBool():
            self.output_model_part.AddNodalSolutionStepVariable(SW.HEIGHT)
            self.output_model_part.AddNodalSolutionStepVariable(KM.MOMENTUM)
            self.output_model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
            self.output_model_part.AddNodalSolutionStepVariable(KM.ACCELERATION)
            self.output_model_part.AddNodalSolutionStepVariable(SW.VERTICAL_VELOCITY)
            self.output_model_part.AddNodalSolutionStepVariable(SW.TOPOGRAPHY)

        # domain_size = self.sw_model_part.ProcessInfo[KM.DOMAIN_SIZE]
        domain_size = 3
        element_name = "Element{}D3N".format(domain_size)
        condition_name = "LineCondition{}D2N".format(domain_size)
        KM.DuplicateMeshModeler(self.sw_model_part).GenerateMesh(
            self.output_model_part, element_name, condition_name)
        self.output_model_part.ProcessInfo[KM.DOMAIN_SIZE] = domain_size

    def _SetOutputProcessInfo(self):
        time = self.sw_model_part.ProcessInfo[KM.TIME]
        step = self.sw_model_part.ProcessInfo[KM.STEP]
        self.output_model_part.ProcessInfo[KM.TIME] = time
        self.output_model_part.ProcessInfo[KM.STEP] = step

    def _MapToOutputModelPart(self):

        if self.settings["store_historical_database"].GetBool():
            for variable in self.variables:

                KM.VariableUtils().CopyModelPartNodalVar(
                    variable,
                    self.sw_model_part,
                    self.output_model_part,
                    0)
        if self.settings["store_non_historical_database"].GetBool():
            for variable in self.variables:
                KM.VariableUtils().CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar(
                    variable, variable,
                    self.sw_model_part,
                    self.output_model_part,
                    KM.Flags(), False)