import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
import KratosMultiphysics.PfemFluidDynamicsApplication as PFEM
import KratosMultiphysics.DelaunayMeshingApplication as Delaunay
from KratosMultiphysics.HDF5Application import read_model_part_from_hdf5_process
from KratosMultiphysics.HDF5Application import single_mesh_temporal_input_process

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DepthIntegrationInputProcess(model, settings["Parameters"])

class DepthIntegrationInputProcess(KM.OutputProcess):
    """DepthIntegrationInputProcess

    Read the depth integrated values from an HDF5 file and set them as boundary conditions.
    """

    @staticmethod
    def GetDefaultParameters():
        return KM.Parameters("""{
            "model_part_name"           : "",
            "interface_model_part_name" : "interface_model_part",
            "read_historical_database"  : false,
            "interval"                  : [0.0,"End"],
            "swap_yz_axis"              : false,
            "ignore_vertical_component" : true,
            "file_settings"             : {}
        }""")

    def __init__(self, model, settings):
        """The constructor of the DepthIntegrationInputProcess"""

        KM.OutputProcess.__init__(self)
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model_part = model[settings["model_part_name"].GetString()]
        self.read_historical = settings["read_historical_database"].GetBool()
        self.interval = KM.IntervalUtility(settings)
        interface_model_part_name = settings["interface_model_part_name"].GetString()
        if model.HasModelPart(interface_model_part_name):
            raise Exception("DepthIntegrationInputProcess: There is an existing interface model part with name '{}'".format(interface_model_part_name))
        self.interface_model_part = model.CreateModelPart(interface_model_part_name)
        self.swap_yz_axis = settings["swap_yz_axis"].GetBool()
        self.ignore_vertical_component = settings["ignore_vertical_component"].GetBool()

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

        self.hdf5_read = read_model_part_from_hdf5_process.Factory(hdf5_process_settings.Clone(), model)
        self.hdf5_process = single_mesh_temporal_input_process.Factory(hdf5_process_settings, model)

    def Check(self):
        '''Check the processes.'''
        self.hdf5_read.Check()
        self.hdf5_process.Check()

    def ExecuteInitialize(self):
        '''Read the interface_model_part and set the variables.'''
        self.hdf5_read.ExecuteInitialize()
        self._CheckInputVariables()
        self._MapToBoundaryCondition()

    def ExecuteInitializeSolutionStep(self):
        '''Set the variables in interface_model_part at the current time step.'''
        self.hdf5_process.ExecuteInitializeSolutionStep() # look for the current time
        self._CheckInputVariables()
        self._MapToBoundaryCondition()

        ###
        print("ExecuteInitializeSolutionStep")
        for node in self.interface_model_part.Nodes:
            print(node.GetValue(KM.VELOCITY))

    def _CheckInputVariables(self):
        if self.swap_yz_axis:
            if self.read_historical:
                SW.ShallowWaterUtilities().SwapYZComponents(KM.MOMENTUM, self.interface_model_part.Nodes)
                SW.ShallowWaterUtilities().SwapYZComponents(KM.VELOCITY, self.interface_model_part.Nodes)
            else:
                SW.ShallowWaterUtilities().SwapYZComponentsNonHistorical(KM.MOMENTUM, self.interface_model_part.Nodes)
                SW.ShallowWaterUtilities().SwapYZComponentsNonHistorical(KM.VELOCITY, self.interface_model_part.Nodes)
        if self.ignore_vertical_component:
            if self.read_historical:
                KM.VariableUtils().SetVariableToZero(KM.MOMENTUM_Z, self.interface_model_part.Nodes)
                KM.VariableUtils().SetVariableToZero(KM.VELOCITY_Z, self.interface_model_part.Nodes)
            else:
                KM.VariableUtils().SetNonHistoricalVariableToZero(KM.MOMENTUM_Z, self.interface_model_part.Nodes)
                KM.VariableUtils().SetNonHistoricalVariableToZero(KM.VELOCITY_Z, self.interface_model_part.Nodes)

    def _MapToBoundaryCondition(self):
        pass
