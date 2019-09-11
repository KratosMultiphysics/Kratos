from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

# Other imports
import numpy as np
import json
import os

#def Create(solver_settings):
#    return SDoFSolver(solver_settings)

class SDoFSolver(object):
    def __init__(self, input_file_name):

        #input_file_name = self.solver_settings["input_file"]
        if not input_file_name.endswith(".json"):
            input_file_name += ".json"

        with open(input_file_name,'r') as ProjectParameters:
            parameters = json.load(ProjectParameters)

        default_settings = {
                "system_parameters":{
                    "mass"      : 100.0,
                    "stiffness" : 4000.0,
                    "damping"   : 6000.0
                },
                "time_integration_parameters":{
                    "alpha_m"   : -0.5,
                    "alpha_f"   : -0.5,
                    "start_time": 0.0,
                    "time_step" : 0.05
                },
                "initial_values":{
                    "displacement"  : 0.0,
                    "velocity"      : 0.0,
                    "acceleration"  : 0.0
                },
                "boundary_conditions":{
                    "external_load" : 0.0
                },
                "solver_parameters": {
                    "buffer_size"   : 3
                },
                "output_parameters":{
                    "file_name" : "sdof_solver/results_sdof.dat"
                }}

        RecursivelyValidateAndAssignDefaults(default_settings, parameters)

        self.mass = parameters["system_parameters"]["mass"]
        self.stiffness = parameters["system_parameters"]["stiffness"]
        self.damping = parameters["system_parameters"]["damping"]

        self.alpha_m = parameters["time_integration_parameters"]["alpha_m"]
        self.alpha_f = parameters["time_integration_parameters"]["alpha_f"]
        self.delta_t = parameters["time_integration_parameters"]["time_step"]
        self.start_time = parameters["time_integration_parameters"]["start_time"]

        self.initial_displacement = parameters["initial_values"]["displacement"]
        self.initial_velocity = parameters["initial_values"]["velocity"]

        self.force = parameters["boundary_conditions"]["external_load"]

        #calculate initial acceleration
        factor = self.force - self.stiffness * self.initial_displacement
        self.initial_acceleration = (1/self.mass) * factor

        beta = 0.25 * (1- self.alpha_m + self.alpha_f)**2
        gamma =  0.50 - self.alpha_m + self.alpha_f

        self.LHS = np.array([[1.0, 0.0, -self.delta_t**2 * beta],
                              [0.0, 1.0, -self.delta_t * gamma],
                              [(1-self.alpha_f)*self.stiffness,
                               (1-self.alpha_f) * self.damping,
                               (1-self.alpha_m) * self.mass]])

        self.RHS_matrix = np.array([[1.0, self.delta_t, self.delta_t**2 * (0.5 - beta)],
                                     [0.0, 1.0, self.delta_t*(1-gamma)],
                                     [-self.alpha_f * self.stiffness,
                                      -self.alpha_f * self.damping,
                                      -self.alpha_m * self.mass]])

        self.buffer_size = parameters["solver_parameters"]["buffer_size"]

        self.output_file_name = parameters["output_parameters"]["file_name"]


    def Initialize(self):
        self.x = np.zeros((3, self.buffer_size))

        initial_values = np.array([self.initial_displacement,
                                   self.initial_velocity,
                                   self.initial_acceleration])
        self.dx = initial_values

        #x and dx contain: [displacement, velocity, acceleration]

        if os.path.isfile(self.output_file_name):
            os.remove(self.output_file_name)

        #apply external load as an initial impulse
        self.load_vector = np.array([0,
                                     0,
                                     self.force])

        self.time = self.start_time

    def OutputSolutionStep(self):
        with open(self.output_file_name, "a") as results_sdof:
            #outputs displacements
            results_sdof.write(str(self.time) + "\t" + str(self.dx[0]) + "\n")

    def AdvanceInTime(self, current_time):
        # similar to the Kratos CloneTimeStep function
        # advances values along the buffer axis (so rolling columns) using numpy's roll
        self.x = np.roll(self.x,1,axis=1)
        # overwriting at the buffer_idx=0 the newest values
        #buffer_idx = 0
        #self.x[:,buffer_idx] = self.dx
        self.x[:,0] = self.dx

        self.time = current_time + self.delta_t
        return self.time

    def SolveSolutionStep(self):
        ## PMT: ToDo: check with Andreas what this intends to do
        #b = self.RHS_matrix @ self.x[:,0]
        b = np.dot(self.RHS_matrix,self.x[:,0])

        #external load only for testing
        b += self.load_vector
        self.dx = np.linalg.solve(self.LHS, b)

    def GetBufferSize(self):
        return self.buffer_size

    def GetDeltaTime(self):
        return self.delta_t

    def GetSolutionStepValue(self, identifier, buffer_idx=0):
        if identifier == "DISPLACEMENT":
            return self.x[:,buffer_idx][0]
        elif identifier == "VELOCITY":
            return self.x[:,buffer_idx][1]
        elif identifier == "ACCELERATION":
            return self.x[:,buffer_idx][2]
        else:
            raise Exception("Identifier is unknown!")

    def SetSolutionStepValue(self, identifier, value, buffer_idx=0):
        if identifier == "DISPLACEMENT":
            self.x[:,buffer_idx][0] = value
        elif identifier == "VELOCITY":
            self.x[:,buffer_idx][1] = value
        elif identifier == "ACCELERATION":
            self.x[:,buffer_idx][2] = value
        else:
            raise Exception("Identifier is unknown!")

    def SetData(self, identifier, data):
        if identifier == "LOAD":
            # last index is the external force
            self.load_vector[-1] = data
        elif identifier == "DISPLACEMENT":
            # first index is displacement
            # maybe use buffer index
            self.SetSolutionStepValue("DISPLACEMENT", data, 0)
        else:
            raise Exception("Identifier is unknown!")

    def GetData(self, identifier):
        if identifier == "LOAD":
            # last index is the external force
            return self.load_vector[-1]
        elif identifier == "DISPLACEMENT":
            # first index is displacement
            return self.GetSolutionStepValue("DISPLACEMENT", 0)
        else:
            raise Exception("Identifier is unknown!")

    def _GetIOName(self):
        return "sdof"

    def _Name(self):
        return self.__class__.__name__

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