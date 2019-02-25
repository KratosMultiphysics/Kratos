from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from time_integration_base_scheme import TimeIntegrationBaseScheme

# Importing tools
from co_simulation_tools import ValidateAndAssignDefaults

# Other imports
import numpy as np
import json
import os

def CreateScheme(scheme_settings):
    return TimeIntegrationGeneralizedAlphaScheme(scheme_settings)

class TimeIntegrationGeneralizedAlphaScheme(TimeIntegrationBaseScheme):
    """
    A single-degree-of-freedom SDoF model

    Using for testing of the MDoF solver
    """
    def __init__(self, scheme_settings):

        default_settings = {
                "type"      : "generalized_alpha",
                "settings"  : {}
            }

        ValidateAndAssignDefaults(default_settings, scheme_settings)

        # validate, assign and remove custom settings
        key = "p_inf"
        if key in scheme_settings["settings"].keys():
            pInf = scheme_settings["settings"]["p_inf"]
            scheme_settings["settings"].pop("p_inf")
        else:
            err_msg  = 'The item with name "' + key
            err_msg += '" is not present in the settings\n'
            raise Exception(err_msg)

        # add buffer size - this is not user-specified
        # each derived scheme specifies it
        scheme_settings["settings"].update({"buffer_size":2})

        # base scheme settings
        super(TimeIntegrationGeneralizedAlphaScheme, self).__init__(scheme_settings["settings"])

        # custom scheme settiungs
        self.alphaM = (2.0 * pInf - 1.0) / (pInf + 1.0)
        self.alphaF = pInf / (pInf + 1.0)
        self.beta = 0.25 * (1 - self.alphaM + self.alphaF)**2
        self.gamma = 0.5 - self.alphaM + self.alphaF

        # coefficients for LHS
        self.a1h = (1.0 - self.alphaM) / (self.beta * self.dt**2)
        self.a2h = (1.0 - self.alphaF) * self.gamma / (self.beta * self.dt)
        self.a3h = 1.0 - self.alphaF

        # coefficients for mass
        self.a1m = self.a1h
        self.a2m = self.a1h * self.dt
        self.a3m = (1.0 - self.alphaM - 2.0 * self.beta) / (2.0 * self.beta)

        #coefficients for damping
        self.a1b = (1.0 - self.alphaF) * self.gamma / (self.beta * self.dt)
        self.a2b = (1.0 - self.alphaF) * self.gamma / self.beta - 1.0
        self.a3b = (1.0 - self.alphaF) * (0.5 * self.gamma / self.beta - 1.0) * self.dt

        # coefficient for stiffness
        self.a1k = -1.0 * self.alphaF

        # coefficients for velocity update
        self.a1v = self.gamma / (self.beta * self.dt)
        self.a2v = 1.0 - self.gamma / self.beta
        self.a3v = (1.0 - self.gamma / (2 * self.beta)) * self.dt

        # coefficients for acceleration update
        self.a1a = self.a1v / (self.dt * self.gamma)
        self.a2a = -1.0 / (self.beta * self.dt)
        self.a3a = 1.0 - 1.0 / (2.0 * self.beta)

    def Initialize(self, model):
        """
        """
        # call function from base
        super(TimeIntegrationGeneralizedAlphaScheme, self).Initialize(model)

		# overwrite with scheme-specific values
        self.buffer[3,0,:] = np.dot(model.m,self.buffer[2,0,:]) + np.dot(model.b,self.buffer[1,0,:]) + np.dot(model.k,self.buffer[0,0,:])
        self.buffer[3,1,:] = np.dot(model.m,self.buffer[2,1,:]) + np.dot(model.b,self.buffer[1,1,:]) + np.dot(model.k,self.buffer[0,1,:])

    def _AssembleLHS(self, model):
        """
        """
        return self.a1h * model.m + self.a2h * model.b + self.a3h * model.k

    def _AssembleRHS(self, model):
        """
        """
        # PMT -> F = ... -> self.f1 is actually for this step a new external force which is read in
        # so use SetSolutionStep() in case there is an external force
        # should be AssembleRHS(self, model, f1)
        #F = (1.0 - self.alphaF) * f1 + self.alphaF * self.f0

        f = (1.0 - self.alphaF) * self.force + self.alphaF * self.buffer[3,1,:]
        RHS = np.dot(model.m,(self.a1m * self.buffer[0,1,:] + self.a2m * self.buffer[1,1,:] + self.a3m * self.buffer[2,1,:]))
        RHS += np.dot(model.b,(self.a1b * self.buffer[0,1,:] + self.a2b * self.buffer[1,1,:] + self.a3b * self.buffer[2,1,:]))
        RHS += np.dot(self.a1k * model.k, self.buffer[0,1,:]) + f

        self.buffer[3,0,:] = self.force
        return RHS

    def UpdateDerivedValues(self):
        """
        """
        self.buffer[1,0,:] = self.a1v * (self.buffer[0,0,:] - self.buffer[0,1,:]) + self.a2v * self.buffer[1,1,:] + self.a3v * self.buffer[2,1,:]
        self.buffer[2,0,:] = self.a1a * (self.buffer[0,0,:] - self.buffer[0,1,:]) + self.a2a * self.buffer[1,1,:] + self.a3a * self.buffer[2,1,:]