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
    def __init__(self, settings, solver_wrappers, process_info, data_communicator):
        super().__init__(settings, process_info, data_communicator)
        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)
        self.variable = self.interface_data.variable
        self.dimension = self.interface_data.dimension
        self.use_transpose = self.settings["use_transpose"].GetBool()

    def Execute(self):
        
        if not self.interface_data.IsDefinedOnThisRank(): return

        process_info = self.interface_data.GetModelPart().ProcessInfo
        time = process_info[KM.TIME]

        if not KM.IntervalUtility(self.settings).IsInInterval(time):
            if self.echo_level > 0:
                cs_tools.cs_print_info("Elemental_data_to_Nodal_data", "Skipped, not in interval")
            return
        
        model_part_interface:KM.ModelPart = self.interface_data.GetModelPart()

        if self.use_transpose:
            ConversionUtilities.ConvertNodalDataToElementalDataTranspose(model_part_interface, self.variable, self.variable)
        else:   
            ConversionUtilities.ConvertNodalDataToElementalDataDirect(model_part_interface, self.variable, self.variable)
        if self.echo_level > 0:
            cs_tools.cs_print_info("ConvertNodalDataToElementalData", "Done")

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"    : "UNSPECIFIED",
            "data_name" : "UNSPECIFIED",
            "use_transpose" : true,
            "interval"  : [0.0, 1e30]
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults


