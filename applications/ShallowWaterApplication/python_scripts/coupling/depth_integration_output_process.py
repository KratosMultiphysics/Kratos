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
            "store_historical_database" : false,
            "direction_of_integration"  : [0.0, 0.0, 1.0],
            "interval"                  : [0.0,"End"],
            "file_settings"             : {},
            "output_time_settings"      : {}
        }""")

    def __init__(self, model, settings):
        """The constructor of the DepthIntegrationOutputProcess"""

        KM.OutputProcess.__init__(self)
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.volume_model_part = model[settings["volume_model_part_name"].GetString()]
        self.interface_model_part = model[settings["interface_model_part_name"].GetString()]
        self.store_historical = settings["store_historical_database"].GetBool()
        self.interval = KM.IntervalUtility(settings)

        integration_settings = KM.Parameters()
        integration_settings.AddValue("volume_model_part_name", settings["volume_model_part_name"])
        integration_settings.AddValue("interface_model_part_name", settings["interface_model_part_name"])
        integration_settings.AddValue("direction_of_integration", settings["direction_of_integration"])
        integration_settings.AddValue("store_historical_database", settings["store_historical_database"])

        self.integration_process = SW.DepthIntegrationProcess(model, integration_settings)

        hdf5_settings = KM.Parameters()
        hdf5_settings.AddValue("model_part_name", settings["interface_model_part_name"])
        hdf5_settings.AddValue("file_settings", settings["file_settings"])
        hdf5_settings.AddValue("output_time_settings", settings["output_time_settings"])
        data_settings = KM.Parameters("""{"list_of_variables" : ["MOMENTUM","VELOCITY","HEIGHT"]}""")
        if self.store_historical:
            hdf5_settings.AddValue("nodal_solution_step_data_settings", data_settings)
        else:
            hdf5_settings.AddValue("nodal_data_value_settings", data_settings)
        hdf5_process_settings = KM.Parameters()
        hdf5_process_settings.AddValue("Parameters", hdf5_settings)
        
        self.hdf5_process = single_mesh_temporal_output_process.Factory(hdf5_process_settings, model)

    def Check(self):
        '''Check the processes.'''
        self.integration_process.Check()
        self.hdf5_process.Check()

    def ExecuteInitialize(self):
        '''Initialize the variables.'''
        if not self.store_historical:
            KM.VariableUtils().SetNonHistoricalVariableToZero(KM.VELOCITY, self.interface_model_part.Nodes)
            KM.VariableUtils().SetNonHistoricalVariableToZero(KM.MOMENTUM, self.interface_model_part.Nodes)
            KM.VariableUtils().SetNonHistoricalVariableToZero(SW.HEIGHT, self.interface_model_part.Nodes)

    def ExecuteBeforeSolutionLoop(self):
        '''Write the interface_model_part in HDF5 format.'''
        self.hdf5_process.ExecuteBeforeSolutionLoop()

    def ExecuteBeforeOutputStep(self):
        '''Perform the depth integration over the interface_model_part.'''
        self.integration_process.Execute()

    def IsOutputStep(self):
        '''IsOutputStep.'''
        # return self.hdf5_process.IsOutputstep()
        return True

    def PrintOutput(self):
        '''Print the integrated variables from interface_model_part.'''
        self.hdf5_process.ExecuteFinalizeSolutionStep()
