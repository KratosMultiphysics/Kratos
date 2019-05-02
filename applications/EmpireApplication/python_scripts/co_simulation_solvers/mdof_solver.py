from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from co_simulation_solvers.co_simulation_base_solver import CoSimulationBaseSolver
from co_simulation_tools import RecursivelyValidateAndAssignDefaults

# Importing models and schemes
from co_simulation_solvers.mdof.solver_models import *
from co_simulation_solvers.mdof.schemes import *

# Other imports
import numpy as np
import json
import os

def CreateSolver(cosim_solver_settings, level):
    return MDoFSolver(cosim_solver_settings, level)

class MDoFSolver(CoSimulationBaseSolver):
    def __init__(self, cosim_solver_settings, level):
        super(MDoFSolver, self).__init__(cosim_solver_settings, level)

        input_file_name = self.cosim_solver_settings["input_file"]
        if not input_file_name.endswith(".json"):
            input_file_name += ".json"

        with open(input_file_name,'r') as ProjectParameters:
            parameters = json.load(ProjectParameters)

        # creating model using a certain module
        model_type = parameters["model_parameters"]["type"]
        model_module = __import__("mdof_" + model_type + "_model")
        self.model = model_module.CreateModel(parameters["model_parameters"])

        # creating model using a certain module
        scheme_type = parameters["time_integration_scheme_parameters"]["type"]
        scheme_module = __import__("time_integration_" + scheme_type + "_scheme")
        # adding the number of dofs for the constructor of the scheme
        nr_of_dofs = len(self.model.u0)
        if "settings" in parameters["time_integration_scheme_parameters"]:
            parameters["time_integration_scheme_parameters"]["settings"].update({"nr_of_dofs":len(self.model.u0)})
        else:
            parameters["time_integration_scheme_parameters"].update({"settings":{"nr_of_dofs":len(self.model.u0)}})
        self.scheme = scheme_module.CreateScheme(parameters["time_integration_scheme_parameters"])

        default_solver_settings = {
                    "buffer_size": 3
                }
        default_output_settings = {
                    "file_name": "results_mdof.dat"
                }

        # solver settings
        RecursivelyValidateAndAssignDefaults(default_solver_settings, parameters["solver_parameters"])
        # output settings
        RecursivelyValidateAndAssignDefaults(default_output_settings, parameters["output_parameters"])

        self.buffer_size = parameters["solver_parameters"]["buffer_size"]
        self.output_file_name = parameters["output_parameters"]["file_name"]

        # PMT this might not be needed
        self.load_vector = None

        # 1st dimension: variables: disp, acc, vel
        # 2nd dimension: buffer size -> specified by user
        # see distinction to the buffer size for the scheme which is scheme-specific
        # 3rd dimension: number of dofs
        self.buffer = np.zeros((3, self.buffer_size, nr_of_dofs))

    def Initialize(self):
        if os.path.isfile(self.output_file_name):
            os.remove(self.output_file_name)

        # initialize scheme parameters
        self.scheme.Initialize(self.model)

        # PMT is this needed?
        self.load_vector = self.scheme.GetLoad()

    def Predict(self):
        return self.scheme.Predict()

    def OutputSolutionStep(self):
        with open(self.output_file_name, "a") as results_sdof:
            # outputs (generic) displacements for each dof
            results_sdof.write(str(self.time) + "\t" + " ".join(str(value) for value in self.buffer[0,0,:]) + "\n")

    def AdvanceInTime(self, current_time):
        # similar to the Kratos CloneTimeStep function
        # advances values along the buffer axis (so rolling columns) using numpy's roll

        self.buffer = np.roll(self.buffer,1,axis=1)
        # overwriting at the buffer_idx=0 the newest values
        self.buffer[0,0,:] = self.scheme.GetDisplacement()
        self.buffer[1,0,:] = self.scheme.GetVelocity()
        self.buffer[2,0,:] = self.scheme.GetAcceleration()

        # update displacement, velocity and acceleration
        self.scheme.AdvanceScheme()
        self.time = current_time + self.scheme.dt

        return self.time

    def SolveSolutionStep(self):
        # sys of eq reads: LHS * u1 = RHS
        self.scheme.Solve(self.model)
        # from u1 determine v1, a1
        self.scheme.UpdateDerivedValues()

    def GetBufferSize(self):
        return self.buffer_size

    def GetDeltaTime(self):
        return self.scheme.dt

    def GetSolutionStepValue(self, identifier, buffer_idx=0):
        # PMT: what should be get: the buffer value in solver or scheme?
        if identifier == "DISPLACEMENT":
            return self.buffer[0,buffer_idx,:]
        elif identifier == "VELOCITY":
            return self.buffer[1,buffer_idx,:]
        elif identifier == "ACCELERATION":
            return self.buffer[2,buffer_idx,:]
        else:
            raise Exception("Identifier is unknown!")

    def SetSolutionStepValue(self, identifier, value, buffer_idx=0):
        # PMT: what should be set: the buffer value in solver or scheme?
        if identifier == "DISPLACEMENT":
            self.buffer[0,buffer_idx,:] = value
        elif identifier == "VELOCITY":
            self.buffer[1,buffer_idx,:] = value
        elif identifier == "ACCELERATION":
            self.buffer[2,buffer_idx,:] = value
        else:
            raise Exception("Identifier is unknown!")

    def SetData(self, identifier, data):
        #
        # PMT: check if this syntax is still ok for MDoF
        #

        if identifier == "LOAD":
            # last index is the external force
            self.load_vector[-1] = data
        elif identifier == "DISPLACEMENT":
            # first index is displacement
            # maybe use buffer index
            self.SetSolutionStepValue("DISPLACEMENT", data)
        else:
            raise Exception("Identifier is unknown!")

    def GetData(self, identifier):
        #
        # PMT: check if this syntax is still ok for MDoF
        #

        if identifier == "LOAD":
            # last index is the external force
            return self.load_vector[-1]
        elif identifier == "DISPLACEMENT":
            # first index is displacement
            return self.GetSolutionStepValue("DISPLACEMENT")
        else:
            raise Exception("Identifier is unknown!")

    def _GetIOName(self):
        return "mdof"

    def _Name(self):
        return self.__class__.__name__