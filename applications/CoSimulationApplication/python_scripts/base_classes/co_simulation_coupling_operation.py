from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

class CoSimulationCouplingOperation(object):
    """Baseclass for the coupling operations used for CoSimulation
    This class can be used to customize the behavior of the CoSimulation,
    by providing a large interface and access to the solvers/models
    """
    def __init__(self, settings):
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self._GetDefaultSettings())
        self.echo_level = self.settings["echo_level"].GetInt()

    def Initialize(self):
        pass

    def Finalize(self):
        pass


    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass


    def InitializeCouplingIteration(self):
        pass

    def FinalizeCouplingIteration(self):
        pass


    def Execute(self):
        raise NotImplementedError('"Execute" is not implemented for {}!'.format(self._ClassName))


    def PrintInfo(self):
        pass

    def Check(self):
        pass

    @classmethod
    def _ClassName(cls):
        return cls.__name__

    @classmethod
    def _GetDefaultSettings(cls):
        return KM.Parameters("""{
            "type"       : "UNSPECIFIED",
            "echo_level" : 0
        }""")
