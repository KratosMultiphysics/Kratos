# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation
from KratosMultiphysics.CoSimulationApplication import ConversionUtilities

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(*args):
    return ElementalToNodalData(*args)

class ElementalToNodalData(CoSimulationCouplingOperation):
    """This operation maps the Elemental Data to Nodal Data for a given ModelPart
    """
    def __init__(self, settings, solver_wrappers, process_info):
        super().__init__(settings, process_info)
        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)

    def Execute(self):
        model_part_interface = self.interface_data.GetModelPart()
        ConversionUtilities.ConvertElementalDataToNodalData(model_part_interface)

        if self.echo_level > 0:
            cs_tools.cs_print_info("Elemental to Nodal Mapping Operation done")

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



