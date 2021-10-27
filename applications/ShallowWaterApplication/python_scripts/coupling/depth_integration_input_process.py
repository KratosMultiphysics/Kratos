import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
from KratosMultiphysics.HDF5Application import single_mesh_temporal_input_process

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
            "model_part_name"           : "",
            "interface_model_part_name" : "interface_model_part",
            "read_historical_database"  : false,
            "interval"                  : [0.0,"End"],
            "invert_yz_axis"            : false,
            "ignore_vertical_component" : true,
            "file_settings"             : {}
        }""")

    def __init__(self, model, settings):
        """The constructor of the DepthIntegrationOutputProcess"""

        KM.OutputProcess.__init__(self)

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model_part = model[self.settings["volume_model_part_name"].GetString()]
        self.read_historical = settings["read_historical_database"].GetBool()
        self.interval = KM.IntervalUtility(settings)
        interface_model_part_name = self.settings["volume_model_part_name"].GetString()
        if model.HasModelPart(interface_model_part_name):
            raise Exception("DepthIntegrationInputProcess: There is an existing interface model part with name '{}'".format(interface_model_part_name))
        self.interface_model_part = model.CreateModelPart(interface_model_part_name)

        hdf5_settings = KM.Parameters()
        hdf5_settings.AddValue("model_part_name", settings["interface_model_part_name"])
        hdf5_settings.AddValue("file_settings", settings["file_settings"])
        data_settings = KM.Parameters("""{"list_of_variables" : ["MOMENTUM","VELOCITY","HEIGHT"]}""")
        if self.read_historical:
            hdf5_settings.AddValue("nodal_solution_step_data_settings", data_settings)
        else:
            hdf5_settings.AddValue("nodal_data_value_settings", data_settings)
        hdf5_process_settings = KM.Parameters()
        hdf5_process_settings.AddValue("Parameters", hdf5_settings)
        
        self.hdf5_process = single_mesh_temporal_input_process.Factory(hdf5_process_settings, model)

    def Check(self):
        self.hdf5_process.Check()

    def ExecuteInitialize(self):
        self.hdf5_process.ExecuteInitialize()
        if not self.read_historical:
            KM.VariableUtils().SetNonHistoricalVariableToZero(KM.VELOCITY, self.interface_model_part.Nodes)
            KM.VariableUtils().SetNonHistoricalVariableToZero(KM.MOMENTUM, self.interface_model_part.Nodes)
            KM.VariableUtils().SetNonHistoricalVariableToZero(SW.HEIGHT, self.interface_model_part.Nodes)

    def ExecuteBeforeSolutionLoop(self):
        self.hdf5_process.ExecuteBeforeSolutionLoop()

    def ExecuteBeforeSolutionStep(self):
        self.hdf5_process.ExecuteBeforeSolutionStep()
        print(self.interface_model_part)
        for node in self.interface_model_part.Nodes:
            print(node.GetSolutionStepValue(KM.VELOCITY))
