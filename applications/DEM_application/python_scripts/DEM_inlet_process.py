from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics.DEMApplication import *


class DEM_Inlet_Process:

    # NOTE: this function is designed with the idea that different derived operations are allowed to have different constructors.
    # nevertheless for ease of usage, we propose a rather wide range of "default parameters" which should be sufficient
    # to implement most user operations
    #
    # model_part -> model_part to which the operation will be applied
    # group_container -> container of the groups. Can be queried for the nodes that correspond to a given group Id
    # group_ids -> python array containing the ids of the groups to which the operation should be applied
    # table_ids --> python array containing the ids of the tables to be applied
    # echo_level -> level of expected echo for the operation: echo_level=0 implies no echo
    def __init__(self, solid_model_part, creator, inlet_model_part, parameters, echo_level=0):
        self.inlet = DEM_Inlet(inlet_model_part)
        self.inlet_model_part = inlet_model_part
        self.solid_model_part = solid_model_part
        self.creator = creator
        self.parameters = parameters
        self.echo_level = echo_level

    def PrintInfo(self):
        print("PrintInfo")
    # this function is designed for being called at the beginning of the computations
    # right after reading the model and the groups

    def ExecuteInitialize(self):
        if(self.echo_level > 0):
            print("Finished ExecuteInitialize for Class", self.PrintInfo())

    # this function is designed for being execute once before the solution loop but after all of the
    # solvers where built
    def ExecuteBeforeSolutionLoop(self):
        if(self.echo_level > 0):
            print("Finished ExecuteBeforeSolutionLoop for Class", self.PrintInfo())

    # this function will be executed at every time step BEFORE performing the solve phase
    def ExecuteInitializeSolutionStep(self):
        if(self.echo_level > 0):
            print("Inserting new DEM elements...")
        self.inlet.CreateElementsFromInletMesh(self.solid_model_part, self.inlet_model_part, self.creator)
        print("Finished ExecuteInitializeSolutionStep for Class", self.PrintInfo())

    # this function will be executed at every time step AFTER performing the solve phase
    def ExecuteFinalizeSolutionStep(self):
        if(self.echo_level > 0):
            print("Finished ExecuteFinalizeSolutionStep for Class", self.PrintInfo())

    # this function will be executed at every time step BEFORE  writing the output
    def ExecuteBeforeOutputStep(self):
        if(self.echo_level > 0):
            print("Finished ExecuteBeforeOutputStep for Class", self.PrintInfo())

    # this function will be executed at every time step AFTER writing the output
    def ExecuteAfterOutputStep(self):
        if(self.echo_level > 0):
            print("Finished ExecuteAfterOutputStep for Class", self.PrintInfo())

    # this function is designed for being called at the beginning of the computations
    # right after reading the model and the groups
    def ExecuteFinalize(self):
        if(self.echo_level > 0):
            print("Finished ExecuteFinalize for Class", self.PrintInfo())
