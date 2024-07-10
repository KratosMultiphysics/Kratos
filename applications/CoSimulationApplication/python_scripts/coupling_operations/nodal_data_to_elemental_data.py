# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation
from KratosMultiphysics.CoSimulationApplication import ConversionUtilities

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(*args):
    return NodalToElementalData(*args)

class NodalToElementalData(CoSimulationCouplingOperation):
    """This operation maps the Nodal Data to Elemental Data for a given ModelPart consistently. The elemental value is an averaged value of Nodal values.
    """
    def __init__(self, settings, solver_wrappers, process_info, data_communicator):
        super().__init__(settings, process_info, data_communicator)
        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)
        self.variable = self.interface_data.variable
        self.dimension = self.interface_data.dimension
        self.consistent = self.settings["consistent"].GetBool()

        if not self.consistent:
            raise RuntimeError("Conservative mapper from nodes on elements is not implemented!")

    def Execute(self):
        if not self.interface_data.IsDefinedOnThisRank(): return

        process_info = self.interface_data.GetModelPart().ProcessInfo
        time = process_info[KM.TIME]

        if not KM.IntervalUtility(self.settings).IsInInterval(time):
            if self.echo_level > 0:
                cs_tools.cs_print_info("Elemental_data_to_Nodal_data", "Skipped, not in interval")
            return

        self.ConvertNodalDataToElementalData()
        if self.echo_level > 0:
            cs_tools.cs_print_info("ConvertNodalDataToElementalData", "Done")

    def ConvertNodalDataToElementalData(self):
        model_part_interface:KM.ModelPart = self.interface_data.GetModelPart()
        for element in model_part_interface.Elements:
            num_of_nodes = len(element.GetNodes())
            if self.dimension == 1:
                value = 0.0
                for node in element.GetNodes():
                    value += node.GetValue(self.variable) / num_of_nodes
            elif self.dimension == 3:
                value = KM.Array3()
                for node in element.GetNodes():
                    value += node.GetValue(self.variable) / num_of_nodes
            element.SetValue(self.variable, value)

    def PrintInfo(self):
        pass

    def Check(self):
        # TODO in case the NORMALS are computed with historical variables then you should check if the var is in the ModelPart
        pass

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"    : "UNSPECIFIED",
            "data_name" : "UNSPECIFIED",
            "consistent": true,
            "interval"  : [0.0, 1e30]
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults



