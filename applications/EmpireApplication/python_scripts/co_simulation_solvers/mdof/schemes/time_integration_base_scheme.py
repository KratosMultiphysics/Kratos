from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing tools
from co_simulation_tools import ValidateAndAssignDefaults

# Other imports
import numpy as np
import json
import os

class TimeIntegrationBaseScheme(object):
    """
    """
    def __init__(self, scheme_settings):

        default_settings = {
                "time_step"     : 0.01,
                "buffer_size"   : 1,
                "nr_of_dofs"    : 2
            }

        ValidateAndAssignDefaults(default_settings, scheme_settings)

        # time step
        self.dt = scheme_settings["time_step"]

        # 1st dimension: variables: disp, acc, vel, force
        # 2nd dimension: buffer size -> needed by the scheme NOT user-specified
        # see distinction to the buffer size for the solver which IS user-specified
        # 3rd dimension: number of dofs
        self.buffer_size = scheme_settings["buffer_size"]
        self.buffer = np.zeros((4,
                                self.buffer_size,
                                scheme_settings["nr_of_dofs"]))

        # u0 -> u(current-)0
        # u1 -> u(current-)1
        # ... i-th...
        # ui -> u(current-)i

        # PMT check if this is needed
        self.force = None

    def Initialize(self, model):
        """
        """
        for idx in range(self.buffer_size):
            self.buffer[0,idx,:] = model.u0
            self.buffer[1,idx,:] = model.v0
            self.buffer[2,idx,:] = model.a0
            self.buffer[3,idx,:] = model.f0

        # PMT check if this is needed
        self.force = model.f0

    def Predict(self):
        """
        """
        return 2.0 * self.buffer[0,0,:] - self.buffer[0,1,:]

    def _AssembleLHS(self, model):
        """
        """
        pass

    def _AssembleRHS(self, model):
        """
        """
        pass

    def Solve(self, model):
        # sys of eq reads: LHS * u0 = RHS
        LHS = self._AssembleLHS(model)
        RHS = self._AssembleRHS(model)
        self.buffer[0,0,:] = np.linalg.solve(LHS, RHS)

    def UpdateDerivedValues(self):
        """
        """
        pass

    def GetDisplacement(self):
        """
        """
        return self.buffer[0,0,:]

    def GetVelocity(self):
        """
        """
        return self.buffer[1,0,:]

    def GetAcceleration(self):
        """
        """
        return self.buffer[2,0,:]

    def GetLoad(self):
        """
        """
        return self.buffer[3,0,:]

    def GetPreviousDisplacement(self):
        """
        """
        return self.buffer[0,1,:]

    def GetPreviousVelocity(self):
        """
        """
        return self.buffer[1,1,:]

    def GetPreviousAcceleration(self):
        """
        """
        return self.buffer[2,1,:]

    def GetPreviousLoad(self):
        """
        """
        return self.buffer[3,1,:]

    def AdvanceScheme(self):
        for idx in range(self.buffer_size-1):
            self.buffer[0,-1-idx,:] = self.buffer[0,-1-(idx+1),:]
            self.buffer[1,-1-idx,:] = self.buffer[1,-1-(idx+1),:]
            self.buffer[2,-1-idx,:] = self.buffer[2,-1-(idx+1),:]
            self.buffer[3,-1-idx,:] = self.buffer[3,-1-(idx+1),:]