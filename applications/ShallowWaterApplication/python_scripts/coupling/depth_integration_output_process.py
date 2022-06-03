import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
from KratosMultiphysics.HDF5Application import single_mesh_temporal_output_process

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DepthIntegrationOutputProcess(model, settings["Parameters"])

class DepthIntegrationOutputProcess(KM.OutputProcess):
    """DepthIntegrationOutputProcess

    This process performs a depth integration from a Navier-Stokes domain to a shallow water domain.
    The depth integration values are stored in the nodes of the shallow water domain and
    printed in HDF5 format.
    """

    @staticmethod
    def GetDefaultParameters():
        return KM.Parameters("""{
            "volume_model_part_name"    : "",
            "interface_model_part_name" : "",
            "output_model_part_name"    : "",
            "store_historical_database" : false,
            "extrapolate_boundaries"    : false,
            "print_velocity_profile"    : false,
            "interval"                  : [0.0,"End"],
            "file_settings"             : {},
            "output_time_settings"      : {}
        }""")

    def __init__(self, model, settings):
        """The constructor of the DepthIntegrationOutputProcess"""

        KM.OutputProcess.__init__(self)
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.volume_model_part = model[self.settings["volume_model_part_name"].GetString()]
        self.interface_model_part = model[self.settings["interface_model_part_name"].GetString()]
        self.output_model_part = model.CreateModelPart(self.settings["output_model_part_name"].GetString())
        self.interval = KM.IntervalUtility(self.settings)
        self.variables = [KM.VELOCITY, KM.MOMENTUM, SW.HEIGHT]

        if self.volume_model_part.ProcessInfo[KM.DOMAIN_SIZE] == 2:
            self.integration_process = SW.DepthIntegrationProcess2D(model, self._CreateIntegrationParameters())
        else:
            self.integration_process = SW.DepthIntegrationProcess3D(model, self._CreateIntegrationParameters())
        self.hdf5_process = single_mesh_temporal_output_process.Factory(self._CreateHDF5Parameters(), model)


    def Check(self):
        '''Check the processes.'''
        self.integration_process.Check()
        self.hdf5_process.Check()


    def ExecuteInitialize(self):
        '''Initialize the output model part.'''
        self._InitializeOutputModelPart()
        self._SetOutputProcessInfo()
        if not self.settings["store_historical_database"].GetBool():
            for var in self.variables:
                KM.VariableUtils().SetNonHistoricalVariableToZero(var, self.interface_model_part.Nodes)
                KM.VariableUtils().SetNonHistoricalVariableToZero(var, self.output_model_part.Nodes)


    def ExecuteBeforeSolutionLoop(self):
        '''Write the interface model part in HDF5 format.'''
        self.integration_process.Execute()
        self._MapToOutputModelPart()
        self.hdf5_process.ExecuteBeforeSolutionLoop()


    def ExecuteInitializeSolutionStep(self):
        '''Synchronize the ProcessInfo of the output and interface model part.'''
        self._SetOutputProcessInfo()


    def IsOutputStep(self):
        '''IsOutputStep.'''
        # return self.hdf5_process.IsOutputstep()
        return True


    def PrintOutput(self):
        '''Perform the depth integration over the interface model part.'''
        self.integration_process.Execute()
        self._MapToOutputModelPart()
        self.hdf5_process.ExecuteFinalizeSolutionStep()


    def _InitializeOutputModelPart(self):
        if self.settings["store_historical_database"].GetBool():
            self.output_model_part.AddNodalSolutionStepVariable(SW.HEIGHT)
            self.output_model_part.AddNodalSolutionStepVariable(KM.MOMENTUM)
            self.output_model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        domain_size = self.volume_model_part.ProcessInfo[KM.DOMAIN_SIZE]
        element_name = "Element{}D2N".format(domain_size)
        condition_name = "LineCondition{}D2N".format(domain_size)
        KM.DuplicateMeshModeler(self.interface_model_part).GenerateMesh(
            self.output_model_part, element_name, condition_name)
        self.output_model_part.ProcessInfo[KM.DOMAIN_SIZE] = domain_size


    def _MapToOutputModelPart(self):
        if self.settings["store_historical_database"].GetBool():
            for variable in self.variables:
                KM.VariableUtils().CopyModelPartNodalVar(
                    variable,
                    self.interface_model_part,
                    self.output_model_part,
                    0)
        else:
            for variable in self.variables:
                KM.VariableUtils().CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar(
                    variable, variable,
                    self.interface_model_part,
                    self.output_model_part,
                    KM.Flags(), False)


    def _SetOutputProcessInfo(self):
        time = self.interface_model_part.ProcessInfo[KM.TIME]
        step = self.interface_model_part.ProcessInfo[KM.STEP]
        self.output_model_part.ProcessInfo[KM.TIME] = time
        self.output_model_part.ProcessInfo[KM.STEP] = step


    def _CreateIntegrationParameters(self):
        integration_settings = KM.Parameters()
        integration_settings.AddValue("volume_model_part_name", self.settings["volume_model_part_name"])
        integration_settings.AddValue("interface_model_part_name", self.settings["interface_model_part_name"])
        integration_settings.AddValue("store_historical_database", self.settings["store_historical_database"])
        integration_settings.AddValue("extrapolate_boundaries", self.settings["extrapolate_boundaries"])
        integration_settings.AddValue("print_velocity_profile", self.settings["print_velocity_profile"])
        return integration_settings


    def _CreateHDF5Parameters(self):
        hdf5_settings = KM.Parameters()
        hdf5_settings.AddValue("model_part_name", self.settings["output_model_part_name"])
        hdf5_settings.AddValue("file_settings", self.settings["file_settings"])
        hdf5_settings.AddValue("output_time_settings", self.settings["output_time_settings"])
        data_settings = KM.Parameters("""{"list_of_variables" : ["MOMENTUM","VELOCITY","HEIGHT"]}""")
        if self.settings["store_historical_database"].GetBool():
            hdf5_settings.AddValue("nodal_solution_step_data_settings", data_settings)
        else:
            hdf5_settings.AddValue("nodal_data_value_settings", data_settings)
        hdf5_process_settings = KM.Parameters()
        hdf5_process_settings.AddValue("Parameters", hdf5_settings)
        return hdf5_process_settings
