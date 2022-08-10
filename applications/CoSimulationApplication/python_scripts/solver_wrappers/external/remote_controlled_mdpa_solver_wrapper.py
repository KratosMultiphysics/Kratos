# CoSimulation imports
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KratosCoSim
import KratosMultiphysics.StructuralMechanicsApplication # needed for some variables

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.external.remote_controlled_solver_wrapper import RemoteControlledSolverWrapper

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities

def Create(settings, model, solver_name):
    return RemoteControlledWithModalPartSolverWrapper(settings, model, solver_name)

class RemoteControlledWithModalPartSolverWrapper(RemoteControlledSolverWrapper):
    """This class is a generic wrapper for connecting external solvers that are being remote controlled
    """
    def __init__(self, settings, model, solver_name):
        super(RemoteControlledSolverWrapper, self).__init__(settings, model, solver_name)

        model_part_file = self.settings["model_part_file"].GetString()
        self.model_part_name = self.settings["model_part_name"].GetString()

        if self.model.HasModelPart(self.model_part_name):
            self.main_model_part = self.model.GetModelPart(self.model_part_name)
        else:
            self.main_model_part = self.model.CreateModelPart(self.model_part_name)

        model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        KM.ModelPartIO(model_part_file).ReadModelPart(self.main_model_part)

    def Initialize(self):
        super().Initialize()
        interface_config = KM.Parameters("""{}""")
        interface_config.AddEmptyValue("model_part_name").SetString(self.settings["interface_submodel_part"].GetString())

        self.ExportCouplingInterface(interface_config)

    def ExportCouplingInterface(self, interface_config):
        self.__SendControlSignal("ImportMesh", interface_config) # TODO this can also be geometry at some point
        interface_json = {"model_part_name" : interface_config["model_part_name"].GetString()}
        super(RemoteControlledSolverWrapper,self).ExportCouplingInterface(interface_json)

    def __SendControlSignal(self, signal, settings=None):
        data_config = {
            "type"           : "control_signal",
            "identifier"     : "Structure",
            "control_signal" : signal,
            "settings"       : settings
        }
        self.ExportData(data_config)

    @classmethod
    def _GetDefaultParameters(cls):
        return KM.Parameters("""{
            "type"                    : "",
            "solver_wrapper_settings" : {},
            "io_settings"             : {},
            "data"                    : {},
            "echo_level"              : 0,
            "model_part_file"         : "",
            "model_part_name"         : "",
            "interface_submodel_part" : ""
        }""")
