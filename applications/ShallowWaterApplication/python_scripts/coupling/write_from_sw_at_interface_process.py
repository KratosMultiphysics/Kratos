import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
from KratosMultiphysics.HDF5Application import single_mesh_temporal_output_process
from KratosMultiphysics.ShallowWaterApplication.coupling import depth_integration_output_process as BaseProcess

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return WriteFromSwAtInterfaceProcess(model, settings["Parameters"])

class WriteFromSwAtInterfaceProcess(BaseProcess.DepthIntegrationOutputProcess):
    """WriteFromSwAtInterfaceProcess

    This process stores the varialbes of a SW simulation into specific nodes, 
    used as interface, and printed in HDF5 format.
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
        """The constructor of the WriteFromSwAtInterfaceProcess"""

        KM.OutputProcess.__init__(self)
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.volume_model_part = model[self.settings["volume_model_part_name"].GetString()]
        self.interface_model_part = model[self.settings["interface_model_part_name"].GetString()]
        self.output_model_part = model.CreateModelPart(self.settings["output_model_part_name"].GetString())
        self.interval = KM.IntervalUtility(self.settings)
        self.variables = [KM.VELOCITY, KM.MOMENTUM, SW.HEIGHT]

        if self.volume_model_part.ProcessInfo[KM.DOMAIN_SIZE] == 2:
            self.integration_process = SW.WriteFromSwAtInterfaceProcess2D(model, self._CreateIntegrationParameters())
        else:
            self.integration_process = SW.WriteFromSwAtInterfaceProcess3D(model, self._CreateIntegrationParameters())
        self.hdf5_process = single_mesh_temporal_output_process.Factory(self._CreateHDF5Parameters(), model)

