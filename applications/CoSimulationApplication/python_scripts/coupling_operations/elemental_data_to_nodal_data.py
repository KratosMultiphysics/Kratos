# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation
from KratosMultiphysics.CoSimulationApplication import ConversionUtilities

import KratosMultiphysics.OptimizationApplication as KratosOA

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(*args):
    return ElementalToNodalData(*args)

class ElementalToNodalData(CoSimulationCouplingOperation):
    """"
    ElementalToNodalData is a class that performs the conversion of elemental data to nodal data in a co-simulation context.

    Attributes:
        settings (KM.Parameters): The settings for the operation.
        solver_wrappers (dict): A dictionary of solver wrappers.
        process_info (KM.ProcessInfo): The process information.
        data_communicator (KM.DataCommunicator): The data communicator.
        interface_data (InterfaceData): The interface data object.
        variable (KM.Variable): The variable to be converted.
        use_transpose (bool): Flag indicating whether to use transpose in the conversion.

    Methods:
        __init__(self, settings, solver_wrappers, process_info, data_communicator):
            Initializes the ElementalToNodalData object with the given settings, solver wrappers, process information, and data communicator.

        Execute(self):

        PrintInfo(self):
            Prints information about the operation.

        Check(self):
            Checks the validity of the operation. (Currently not implemented)

        _GetDefaultParameters(cls):
            Returns the default parameters for the operation.
    """

    def __init__(self, settings, solver_wrappers, process_info, data_communicator):
        super().__init__(settings, process_info, data_communicator)
        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)
        self.variable = self.interface_data.variable
        self.use_transpose = self.settings["use_transpose"].GetBool()


    def Execute(self):
        """
        Executes the conversion of elemental data to nodal data.

        This function performs the conversion of elemental data to nodal data
        for the given interface data. It checks if the data is defined on the
        current rank and if the current time is within the specified interval.
        Depending on the `use_transpose` flag, it either uses the transpose (nodal value = elem value / num_nodes or
        direct (simple distribution of the elemental values to nodes) method for the conversion.

        Returns:
            None
        """
        if not self.interface_data.IsDefinedOnThisRank(): return

        process_info = self.interface_data.GetModelPart().ProcessInfo
        time = process_info[KM.TIME]

        if not KM.IntervalUtility(self.settings).IsInInterval(time):
            if self.echo_level > 0:
                cs_tools.cs_print_info("Elemental_data_to_Nodal_data", "Skipped, not in interval")
            return

        model_part_interface = self.interface_data.GetModelPart()
        if self.use_transpose:
            ConversionUtilities.ConvertElementalDataToNodalDataTranspose(model_part_interface, self.variable, self.variable)
        else:
            ConversionUtilities.ConvertElementalDataToNodalDataDirect(model_part_interface, self.variable, self.variable)

        if self.echo_level > 0:
            cs_tools.cs_print_info("Elemental_data_to_Nodal_data", "Done")

    def PrintInfo(self):
        pass

    def Check(self):
        # TODO in case the NORMALS are computed with historical variables then you should check if the var is in the ModelPart
        pass

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"           : "UNSPECIFIED",
            "data_name"        : "UNSPECIFIED",
            "use_transpose"    : true,
            "interval"         : [0.0, 1e30]
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults