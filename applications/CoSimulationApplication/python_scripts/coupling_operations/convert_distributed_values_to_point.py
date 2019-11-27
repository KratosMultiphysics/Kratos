from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

class ConvertDistributedValuesToPoint(CoSimulationCouplingOperation):

    def __init__(self, settings, solver_wrappers):
        super(ComputeNormalsOperation, self).__init__(settings)
        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)

        if self.interface_data.is_scalar_variable:
            self.intermediate_variable = KM.KratosGlobals.GetVariable(self.settings["intermediate_variable_scalar"].GetString())
        else:
            self.intermediate_variable = KM.KratosGlobals.GetVariable(self.settings["intermediate_variable_scalar"].GetString())

    def Execute(self):
        KM.VariableRedistributionUtility.ConvertDistributedValuesToPoint(
            self.interface_data.GetModelPart(),
            self.interface_data.variable,
            self.intermediate_variable)

    def Check(self):
        # currently the redistribution-utility only works with conditions
        num_local_conds = self.interface_data.GetModelPart().GetCommunicator().LocalMesh().NumberOfConditions()
        data_comm = self.interface_data.GetModelPart().GetCommunicator().GetDataCommunicator()
        num_local_conds = data_comm.SumAll(num_local_conds)
        if num_local_conds < 1:
            raise Exception('No conditions found in ModelPart "{}" of interface data "{}" of solver "{}"!'.format(var.Name(), self.interface_data.model_part_name, self.interface_data.name, self.interface_data.solver_name))

        # currently the redistribution-utility only supports historical values on the nodes
        if not self.interface_data.GetModelPart().HasNodalSolutionStepVariable(self.intermediate_variable):
            raise Exception('Variable "{}" is missing as SolutionStepVariable in ModelPart "{}" of interface data "{}" of solver "{}"!'.format(var.Name(), self.interface_data.model_part_name, self.interface_data.name, self.interface_data.solver_name))

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "solver"                       : "UNSPECIFIED",
            "data_name"                    : "UNSPECIFIED",
            "intermediate_variable_scalar" : "NODAL_PAUX",
            "intermediate_variable_vector" : "NODAL_VAUX"
        }""")
        this_defaults.AddMissingParameters(super(DistributePointValuesOperation, cls)._GetDefaultSettings())
        return this_defaults

