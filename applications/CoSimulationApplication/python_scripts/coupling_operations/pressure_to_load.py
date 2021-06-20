# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation
from KratosMultiphysics.CoSimulationApplication import ConversionUtilities

def Create(*args):
    return PressureToLoad(*args)

class PressureToLoad(CoSimulationCouplingOperation):
    """This operation computes the Nodal Load using the Elemental Presure loads in a given ModelPart
    """
    def __init__(self, settings, solver_wrappers, process_info):
        super().__init__(settings, process_info)
        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)

    def Execute(self):
        model_part_interface = self.interface_data.GetModelPart()
        CoSimulationCouplingOperation.ConversionUtilities.ConvertPressureToForces(model_part_interface)

    def PrintInfo(self):
        pass

    def Check(self):
        # TODO in case the NORMALS are computed with historical variables then you should check if the var is in the ModelPart
        pass

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"    : "UNSPECIFIED",
            "data_name" : "UNSPECIFIED"
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults



