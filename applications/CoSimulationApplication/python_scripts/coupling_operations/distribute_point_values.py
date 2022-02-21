# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

def Create(*args):
    return DistributePointValuesOperation(*args)

class DistributePointValuesOperation(CoSimulationCouplingOperation):
    """This operation converts concentrated nodal values (such as nodal loads) into distributed quantities (such as tractions)
    It is the inverse operation of the "ConvertDistributedValuesToPoint" operation
    """

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
        # TODO refactor the utility to allow mixed locations!
        fct_ptr = KM.VariableRedistributionUtility.DistributePointValues
        if self.data_point.location == "node_non_historical":
            fct_ptr = KM.VariableRedistributionUtility.DistributePointValuesNonHistorical

        fct_ptr(
            self.data_point.GetModelPart(),
            self.entities,
            self.data_point.variable,
            self.data_dist.variable,
            self.redistribution_tolerance,
            self.redistribution_iterations)

    def Check(self):
        if self.data_point.model_part_name != self.data_dist.model_part_name:
            raise Exception('The ModelParts must be the same!\n    ModelPart of point-data:       "{}"\n    ModelPart of distributed-data: "{}"'.format(self.data_point.model_part_name, self.data_dist.model_part_name))

        num_entities = len(self.entities)
        data_comm = self.data_point.GetModelPart().GetCommunicator().GetDataCommunicator()
        num_entities = data_comm.SumAll(num_entities)
        if num_entities < 1:
            raise Exception('No entities ("{}") found in ModelPart "{}" of interface data "{}" of solver "{}"!'.format(self.settings["entities"].GetString(), self.data_point.model_part_name, self.data_point.name, self.data_point.solver_name))

        if "node" not in self.data_point.location:
            raise Exception('Only nodal values are supported (and not "{}")!\n    ModelPart "{}" of interface data "{}" of solver "{}"!'.format(self.data_point.location, self.data_point.model_part_name, self.data_point.name, self.data_point.solver_name))

        if self.data_point.location != self.data_dist.location:
            raise Exception('The location of the data must be the same!\n    Location of data in point-data:       "{}"\n    Location of data in distributed-data: "{}"'.format(self.data_point.location, self.data_dist.location))

        if self.data_point.variable_type not in ["Double", "Array"]:
            raise Exception('Only variables of type "Double" or "Array" are supported (and not "{}")!\n    ModelPart "{}" of interface data "{}" of solver "{}"!'.format(self.data_point.variable_type, self.data_point.model_part_name, self.data_point.name, self.data_point.solver_name))

        if self.data_point.variable_type != self.data_dist.variable_type:
            raise Exception('The variable types of the data must be the same!\n    Variable type of data in point-data:       "{}"\n    Variable type of data in distributed-data: "{}"'.format(self.data_point.variable_type, self.data_dist.variable_type))

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
