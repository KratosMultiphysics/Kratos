from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7
import KratosMultiphysics
import co_simulation_ios.co_simulation_io_factory as io_factory
# Importing the base class
from co_simulation_solvers.co_simulation_base_solver import CoSimulationBaseSolver
from co_simulation_tools import RecursivelyValidateAndAssignDefaults

# Other imports
import numpy as np
import json
import os

def CreateSolver(cosim_solver_settings, level):
    return SDoFStaticSolver(cosim_solver_settings, level)

class SDoFStaticSolver(CoSimulationBaseSolver):
    def __init__(self, cosim_solver_settings, level):

        super(SDoFStaticSolver, self).__init__(cosim_solver_settings, level)
        input_file_name = self.cosim_solver_settings["input_file"]
        if not input_file_name.endswith(".json"):
            input_file_name += ".json"

        with open(input_file_name,'r') as ProjectParameters:
            parameters = json.load(ProjectParameters)

        default_settings = {
                "system_parameters":{
                    "stiffness" : 4000.0
                },
                "initial_values":{
                    "displacement"  : 0.0,
                },
                "boundary_conditions":{
                    "external_load" : 5000.0
                },
                "solver_parameters": {
                    "buffer_size"   : 1
                },
                "output_parameters":{
                    "file_name" : "sdof_static_solver/results_sdof.dat"
                }}

        RecursivelyValidateAndAssignDefaults(default_settings, parameters)

        self.stiffness = parameters["system_parameters"]["stiffness"]

        self.initial_displacement = parameters["initial_values"]["displacement"]

        self.force = parameters["boundary_conditions"]["external_load"]

        self.buffer_size = parameters["solver_parameters"]["buffer_size"]

        self.output_file_name = parameters["output_parameters"]["file_name"]

    def Initialize(self):
        self.x = np.zeros((1, self.buffer_size))

        initial_values = self.initial_displacement
        self.dx = initial_values

        #x contain: displacement and dx is not used here

        if os.path.isfile(self.output_file_name):
            os.remove(self.output_file_name)
        self.force

    def OutputSolutionStep(self):
        with open(self.output_file_name, "a") as results_sdof_static:
            #outputs displacements
            # results_sdof_static.write("displacement values are \t" + str(self.x + "\n")
            return
    def AdvanceInTime(self, current_time):
        pass

    def SolveSolutionStep(self):

        b = self.force
        self.x = b/self.stiffness

        print('force = ', b)
        print('displacement = ', self.x)

    def GetBufferSize(self):
        return self.buffer_size

    def GetSolutionStepValue(self, identifier, buffer_idx=0):
        if identifier == "DISPLACEMENT":
            return self.x
        else:
            raise Exception("Identifier is unknown!")

    def SetSolutionStepValue(self, identifier, value, buffer_idx=0):
        if identifier == "DISPLACEMENT":
            self.x= value
        else:
            raise Exception("Identifier is unknown!")

    def GetDataDefinition(self, data_name):
        return self.cosim_solver_settings["data"][data_name]

    def SetEchoLevel(self, level):
        self.echo_level = level

    def SetData(self, identifier, data):
        if identifier == "LOAD":
            self.force = data
        elif identifier == "DISPLACEMENT":
            self.SetSolutionStepValue("DISPLACEMENT", data,0)
        else:
            raise Exception("Identifier is unknown!")

    def GetData(self, identifier):
        if identifier == "LOAD":
            return self.force
        elif identifier == "DISPLACEMENT":
            return self.GetSolutionStepValue("DISPLACEMENT",0)
        else:
            raise Exception("Identifier is unknown!")

    def _GetIOName(self):
        return "sdof"

    def _Name(self):
        return self.__class__.__name__