from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

## This proces sets the value of a vector variable component-by-component.
## In this case, the fixicity is given set by deffault to true.
import sys
from KratosMultiphysics.PfemFluidDynamicsApplication.assign_vector_components_to_nodes_process import AssignVectorComponentsToNodesProcess

def Factory(custom_settings, Model):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignInletProcess(Model, custom_settings["Parameters"])

## All the processes python should be derived from "Process"

class AssignInletProcess(AssignVectorComponentsToNodesProcess):
    def __init__(self, Model, custom_settings ):

        AssignVectorComponentsToNodesProcess.__init__(self, Model, custom_settings)

        self.model_part = Model[custom_settings["model_part_name"].GetString()]

        self.echo_level=0
        self.model_inlet =  KratosPfemFluid.SetInlet(self.model_part,self.echo_level)
        self.model_inlet.Execute()


