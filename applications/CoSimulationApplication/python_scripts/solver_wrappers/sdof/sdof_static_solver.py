# Importing the base class
import KratosMultiphysics
from KratosMultiphysics.CoSimulationApplication.function_callback_utility import GenericCallFunction

# Other imports
import json
import os

class SDoFStaticSolver(object):
    def __init__(self, input_name):

        # mimicing two constructors
        if isinstance(input_name, dict):
            parameters = input_name

        elif isinstance(input_name, str):
            if not input_name.endswith(".json"):
                input_name += ".json"

            with open(input_name,'r') as ProjectParameters:
                parameters = json.load(ProjectParameters)

        else:
            raise Exception("The input has to be provided as a dict or a string")

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
                    "write_output_file": False,
                    "file_name" : "sdof_static_solver/results_sdof.dat"
                }}

        RecursivelyValidateAndAssignDefaults(default_settings, parameters)

        self.stiffness = parameters["system_parameters"]["stiffness"]

        self.initial_displacement = parameters["initial_values"]["displacement"]

        self.force = parameters["boundary_conditions"]["external_load"]

        self.buffer_size = parameters["solver_parameters"]["buffer_size"]

        self.output_file_name = parameters["output_parameters"]["file_name"]

        self.write_output_file = parameters["output_parameters"]["write_output_file"]

    def Initialize(self):
        initial_values = self.initial_displacement
        self.dx = initial_values

        if self.write_output_file:
            if os.path.isfile(self.output_file_name):
                os.remove(self.output_file_name)
            self.InitializeOutput()
        self.time = 0.0

    def InitializeOutput(self):
        with open(self.output_file_name, "w") as results_sdof_static:
            results_sdof_static.write("displacement" + "\n")
        self.OutputSolutionStep()

    def OutputSolutionStep(self):
        if self.write_output_file:
            with open(self.output_file_name, "a") as results_sdof_static:
                #outputs results
                results_sdof_static.write(str(self.dx) + "\n")

    def AdvanceInTime(self, current_time):
        self.time = 0.0
        return self.time

    def SolveSolutionStep(self):
        self.dx = self.force/self.stiffness
        KratosMultiphysics.Logger.PrintInfo('SDoFStaticSolver', 'Force Imported = ', self.force)
        KratosMultiphysics.Logger.PrintInfo('SDoFStaticSolver', 'Structure Stiffness = ', self.stiffness)
        KratosMultiphysics.Logger.PrintInfo('SDoFStaticSolver', 'New Displacement = ', self.dx)

    def CalculateReaction(self, buffer_idx=0):
        reaction = self.stiffness * (self.dx)
        return reaction

    def GetSolutionStepValue(self, identifier, buffer_idx=0):
        if identifier == "DISPLACEMENT":
            return self.dx
        elif identifier == "REACTION":
            return self.CalculateReaction()
        else:
            raise Exception("Identifier is unknown!")

    def SetSolutionStepValue(self, identifier, value, buffer_idx=0):
        if identifier == "DISPLACEMENT":
            self.dx= value
        elif identifier == "LOAD":
            self.force = 0.0
            self.force = value
        elif identifier == "ROOT_POINT_DISPLACEMENT":
            self.root_point_displacement = 0.0
            self.root_point_displacement = value
        else:
            raise Exception("Identifier is unknown!")

def ValidateAndAssignDefaults(defaults, settings, recursive=False):
    for key, val in settings.items():
        # check if the current entry also exists in the defaults
        if not key in defaults.keys():
            err_msg  = 'The item with name "' + key + '" is present in this '
            err_msg += 'settings\nbut NOT in the defaults!\n'
            err_msg += 'settings are:\n'
            err_msg += json.dumps(settings, indent=4)
            err_msg += '\ndefaults are:\n'
            err_msg += json.dumps(defaults, indent=4)
            raise Exception(err_msg)

        # check if the type is the same in the defaults
        if type(settings[key]) != type(defaults[key]):
            err_msg  = 'The type of the item with name "' + key + '" (type: "'
            err_msg += str(type(settings[key]).__name__)+'") in this '
            err_msg += 'settings\nis NOT the same as in the defaults (type: "'
            err_msg += str(type(defaults[key]).__name__)+'")!\n'
            err_msg += 'settings are:\n'
            err_msg += json.dumps(settings, indent=4)
            err_msg += '\ndefaults are:\n'
            err_msg += json.dumps(defaults, indent=4)
            raise Exception(err_msg)

    # loop the defaults and add the missing entries
    for key_d, val_d in defaults.items():
        if key_d not in settings: # add the default in case the setting is not present
            settings[key_d] = val_d
        elif recursive and type(val_d) is dict:
            RecursivelyValidateAndAssignDefaults(val_d, settings[key_d])

def RecursivelyValidateAndAssignDefaults(defaults, settings):
    ValidateAndAssignDefaults(defaults, settings, recursive=True)
