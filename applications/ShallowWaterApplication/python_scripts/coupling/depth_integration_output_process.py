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
            "output_time_settings"      : {},
        }""")

    def __init__(self, model, settings):
        """The constructor of the DepthIntegrationOutputProcess"""

        KM.OutputProcess.__init__(self)

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.volume_model_part = model[self.settings["volume_model_part_name"].GetString()]
        self.interface_model_part = model[self.settings["interface_model_part_name"].GetString()]
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
        data_settings = KM.Parameters("""{"list_of_variables" : ["VELOCITY","HEIGHT"]}""")
        if self.store_historical:
            hdf5_settings.AddValue("nodal_solution_step_data_settings", data_settings)
        else:
            hdf5_settings.AddValue("nodal_data_value_settings", data_settings)
        
        self.hdf5_process = single_mesh_temporal_output_process.Factory(hdf5_settings, model)

    def Check(self):
        self.integration_process.Check()
        self.hdf5_process.Check()

    def ExecuteInitialize(self):
        self.integration_process.ExecuteInitialize()
        self.hdf5_process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        self.integration_process.ExecuteBeforeSolutionLoop()
        self.hdf5_process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.integration_process.ExecuteInitialize()
        self.hdf5_process.ExecuteInitialize()