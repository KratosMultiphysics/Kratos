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
    '''
    ElementalToNodalData is a CoSimulationCouplingOperation that maps Elemental Data to Nodal Data for a given ModelPart.
        Attributes:
        settings (KM.Parameters): Configuration settings for the operation.
        solver_wrappers (dict): Dictionary containing solver wrappers.
        process_info (KM.ProcessInfo): Process information.
        data_communicator (KM.DataCommunicator): Data communicator for parallel execution.
        interface_data (InterfaceData): Interface data object.
        variable (KM.Variable): Variable to be mapped.

    Methods:
        __init__(settings, solver_wrappers, process_info, data_communicator):
            Initializes the ElementalToNodalData operation with the given settings, solver wrappers, process info, and data communicator.
        
        Execute():
            Executes the operation to map Elemental Data to Nodal Data if the current time is within the specified interval.
        
        PrintInfo():
            Prints information about the operation.
        
        Check():
            Checks the validity of the operation. (Currently not implemented)
        
        _GetDefaultParameters():
            Returns the default parameters for the operation.
    '''
    def __init__(self, settings, solver_wrappers, process_info, data_communicator):
        """
        Initializes the ElementalDataToNodalData object.
        Args:
            settings (KratosMultiphysics.Parameters): Configuration settings for the operation.
            solver_wrappers (dict): Dictionary containing solver wrappers.
            process_info (KratosMultiphysics.ProcessInfo): Process information.
            data_communicator (KratosMultiphysics.DataCommunicator): Data communicator for parallel operations.
        """

        super().__init__(settings, process_info, data_communicator)
        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)
        self.variable = self.interface_data.variable


    def Execute(self):
        """
        Executes the conversion of elemental data to nodal data.

        This method performs the following steps:
        1. Checks if the interface data is defined on the current rank. If not, it returns immediately.
        2. Retrieves the current time from the process information of the model part.
        3. Checks if the current time is within the specified interval. If not, it logs a message and returns.
        4. Converts the elemental data to nodal data for the specified variable.
        5. Logs a completion message if the echo level is greater than 0.

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
        ConversionUtilities.ConvertElementalDataToNodalData(model_part_interface, self.variable, self.variable) 

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
            "solver"    : "UNSPECIFIED",
            "data_name" : "UNSPECIFIED",
            "consistent": false,
            "interval"  : [0.0, 1e30]
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults