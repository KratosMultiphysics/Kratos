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
        self.variables = [KM.VELOCITY, KM.MOMENTUM, SW.HEIGHT, SW.VERTICAL_VELOCITY, SW.TOPOGRAPHY]

        if self.volume_model_part.ProcessInfo[KM.DOMAIN_SIZE] == 2:
            self.integration_process = SW.WriteFromSwAtInterfaceProcess2D(model, self._CreateIntegrationParameters())
        else:
            self.integration_process = SW.WriteFromSwAtInterfaceProcess3D(model, self._CreateIntegrationParameters())
        self.hdf5_process = single_mesh_temporal_output_process.Factory(self._CreateHDF5Parameters(), model)

    def _InitializeOutputModelPart(self):
        if self.settings["store_historical_database"].GetBool():
            self.output_model_part.AddNodalSolutionStepVariable(SW.HEIGHT)
            self.output_model_part.AddNodalSolutionStepVariable(KM.MOMENTUM)
            self.output_model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
            self.output_model_part.AddNodalSolutionStepVariable(SW.VERTICAL_VELOCITY)
            self.output_model_part.AddNodalSolutionStepVariable(SW.TOPOGRAPHY)
        domain_size = self.volume_model_part.ProcessInfo[KM.DOMAIN_SIZE]
        element_name = "Element{}D2N".format(domain_size)
        condition_name = "LineCondition{}D2N".format(domain_size)
        KM.DuplicateMeshModeler(self.interface_model_part).GenerateMesh(
            self.output_model_part, element_name, condition_name)
        self.output_model_part.ProcessInfo[KM.DOMAIN_SIZE] = domain_size

    def _CreateHDF5Parameters(self):
        hdf5_settings = KM.Parameters()
        hdf5_settings.AddValue("model_part_name", self.settings["output_model_part_name"])
        hdf5_settings.AddValue("file_settings", self.settings["file_settings"])
        hdf5_settings.AddValue("output_time_settings", self.settings["output_time_settings"])
        data_settings = KM.Parameters("""{"list_of_variables" : ["MOMENTUM","VELOCITY","HEIGHT", "VERTICAL_VELOCITY","TOPOGRAPHY"]}""")
        if self.settings["store_historical_database"].GetBool():
            hdf5_settings.AddValue("nodal_solution_step_data_settings", data_settings)
        else:
            hdf5_settings.AddValue("nodal_data_value_settings", data_settings)
        hdf5_process_settings = KM.Parameters()
        hdf5_process_settings.AddValue("Parameters", hdf5_settings)
        return hdf5_process_settings