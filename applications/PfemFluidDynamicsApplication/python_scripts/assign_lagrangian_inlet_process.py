#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

## This process sets the value of a vector variable component-by-component.
## In this case, the fixicity is given set by default to true.
import sys
from KratosMultiphysics.PfemFluidDynamicsApplication.assign_vector_components_to_nodes_process import AssignVectorComponentsToNodesProcess

def Factory(custom_settings, Model):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignLagrangianInletProcess(Model, custom_settings["Parameters"])

## All the processes python should be derived from "Process"

class AssignLagrangianInletProcess(AssignVectorComponentsToNodesProcess):
    def __init__(self, Model, custom_settings ):

        AssignVectorComponentsToNodesProcess.__init__(self, Model, custom_settings)

        self.model_part = Model[custom_settings["model_part_name"].GetString()]

        self.echo_level=0
        self.model_inlet =  KratosPfemFluid.SetLagrangianInlet(self.model_part,self.echo_level)
        self.model_inlet.Execute()


