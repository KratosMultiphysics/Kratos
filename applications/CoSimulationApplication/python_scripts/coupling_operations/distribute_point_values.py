# Importing the Kratos Library
import KratosMultiphysics as KM
# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

def Create(*args):
    return DistributePointValuesOperation(*args)

class DistributePointValuesOperation(CoSimulationCouplingOperation):

    def __init__(self, settings, solver_wrappers, process_info, data_communicator):
        super().__init__(settings, process_info, data_communicator)

        solver_name = self.settings["solver"].GetString()
        self.data_point = solver_wrappers[solver_name].GetInterfaceData(self.settings["data_point_values"].GetString())
        self.data_dist  = solver_wrappers[solver_name].GetInterfaceData(self.settings["data_distributed_values"].GetString())

        self.redistribution_iterations = self.settings["redistribution_iterations"].GetInt()
        self.redistribution_tolerance  = self.settings["redistribution_tolerance"].GetDouble()

        entities_to_use = self.settings["entities"].GetString()
        if entities_to_use == "conditions":
            self.entities = self.data_point.GetModelPart().Conditions
        elif entities_to_use == "elements":
            self.entities = self.data_point.GetModelPart().Elements
        else:
            raise Exception('"entities" can only be "conditions" or "elements"!')

    def Execute(self):
        # TODO adapt call if non-hist!
        KM.VariableRedistributionUtility.DistributePointValues(
            self.data_point.GetModelPart(),
            self.entities,
            self.data_point.variable,
            self.data_dist.variable,
            self.redistribution_tolerance,
            self.redistribution_iterations)

    def Check(self):
        return
        # check input MPs are the same!
        # Check if have entities
        # Check if is compatible type => nodal or nodal-non-hist
        # Check vars have same type => either both are double or both are array3

        # currently the redistribution-utility only works with conditions
        num_local_conds = self.interface_data.GetModelPart().GetCommunicator().LocalMesh().NumberOfConditions()
        data_comm = self.interface_data.GetModelPart().GetCommunicator().GetDataCommunicator()
        num_local_conds = data_comm.SumAll(num_local_conds)
        if num_local_conds < 1:
            raise Exception('No conditions found in ModelPart "{}" of interface data "{}" of solver "{}"!'.format(var.Name(), self.interface_data.model_part_name, self.interface_data.name, self.interface_data.solver_name))

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"                    : "UNSPECIFIED",
            "data_point_values"         : "UNSPECIFIED",
            "data_distributed_values"   : "UNSPECIFIED",
            "entities"                  : "conditions",
            "redistribution_tolerance"  : 1e-7,
            "redistribution_iterations" : 100
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults

