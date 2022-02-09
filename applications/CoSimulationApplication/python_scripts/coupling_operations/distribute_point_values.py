# Importing the Kratos Library
import KratosMultiphysics as KM
# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

class DistributePointValuesOperation(CoSimulationCouplingOperation):

    def __init__(self, settings, solver_wrappers):
        super(ComputeNormalsOperation, self).__init__(settings)
        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)

        self.redistribution_iterations = self.settings["redistribution_tolerance"].GetDouble()
        self.redistribution_tolerance  = self.settings["redistribution_iterations"].GetInt()

        if self.interface_data.is_scalar_variable:
            self.intermediate_variable = KM.KratosGlobals.GetVariable(self.settings["intermediate_variable_scalar"].GetString())
        else:
            self.intermediate_variable = KM.KratosGlobals.GetVariable(self.settings["intermediate_variable_scalar"].GetString())

    def Execute(self):
        KM.VariableRedistributionUtility.DistributePointValues(
            self.interface_data.GetModelPart(),
            self.intermediate_variable,
            self.interface_data.variable,
            self.redistribution_tolerance,
            self.redistribution_iterations)

    def Check(self):
        # currently the redistribution-utility only works with conditions
        num_local_conds = self.interface_data.GetModelPart().GetCommunicator().LocalMesh().NumberOfConditions()
        data_comm = self.interface_data.GetModelPart().GetCommunicator().GetDataCommunicator()
        num_local_conds = data_comm.SumAll(num_local_conds)
        if num_local_conds < 1:
            raise Exception('No conditions found in ModelPart "{}" of interface data "{}" of solver "{}"!'.format(var.Name(), self.interface_data.model_part_name, self.interface_data.name, self.interface_data.solver_name))

        # currently the redistribution-utility only supports historical values on the nodes
        vars_to_check = [self.intermediate_variable]

        if self.interface_data.is_scalar_variable:
            vars_to_check.append(KM.MAPPER_SCALAR_PROJECTION_RHS)
        else:
            vars_to_check.append(KM.MAPPER_VECTOR_PROJECTION_RHS)

        for var in vars_to_check:
            if not self.interface_data.GetModelPart().HasNodalSolutionStepVariable(var):
                raise Exception('Variable "{}" is missing as SolutionStepVariable in ModelPart "{}" of interface data "{}" of solver "{}"!'.format(var.Name(), self.interface_data.model_part_name, self.interface_data.name, self.interface_data.solver_name))

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "solver"                       : "UNSPECIFIED",
            "data_name"                    : "UNSPECIFIED",
            "redistribution_tolerance"     : 1e-7,
            "redistribution_iterations"    : 100,
            "intermediate_variable_scalar" : "NODAL_PAUX",
            "intermediate_variable_vector" : "NODAL_VAUX"
        }""")
        this_defaults.AddMissingParameters(super(DistributePointValuesOperation, cls)._GetDefaultSettings())
        return this_defaults

