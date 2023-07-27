import math
import KratosMultiphysics
from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyOutletTwoFluidProcess(Model, settings["Parameters"])


class ApplyOutletTwoFluidProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"    : "",
            "penalty_factor": 1.0
        }
        """)
        settings.ValidateAndAssignDefaults(default_settings)

        # Set the OUTLET flag in the outlet model part nodes and conditions
        self.outlet_model_part = Model[settings["model_part_name"].GetString()]
        # self.outlet_model_part.ProcessInfo[KratosFluid.OUTLET_INFLOW_CONTRIBUTION_SWITCH] = False
        domain_size = self.outlet_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        for node in self.outlet_model_part.Nodes:
            node.Set(KratosMultiphysics.OUTLET, True)
        for condition in self.outlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.OUTLET, True)
        root_model_part = self.outlet_model_part.GetRootModelPart()
        for element in root_model_part.Elements:
            for node in element.GetGeometry():
                if node.Is(KratosMultiphysics.OUTLET):
                    element.Set(KratosMultiphysics.OUTLET,True)
                    break




        # Calculate Normal
        self.__CalculateOutletNormal(self.outlet_model_part, domain_size)
        for element in self.outlet_model_part.Elements:
            element.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)

    def __CalculateOutletNormal(self, model_part_name, domain_size):
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplexNonHistorical(
            model_part_name,
            domain_size,
            KratosMultiphysics.DISPLACEMENT)








