from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

class CoSimulationDataTransferOperator(object):
    """Baseclass for the data transfer operators used for CoSimulation
    It transfers data from one interface to another. This can e.g. be mapping or a copy of values.
    """
    def __init__(self, settings):
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self._GetDefaultSettings())
        self.echo_level = self.settings["echo_level"].GetInt()

    def TransferData(self, from_solver_data, to_solver_data, transfer_options):
        raise NotImplementedError("This function has to be implemented in the derived class!")

    def PrintInfo(self):
        pass

    def Check(self):
        pass

    @classmethod
    def _ClassName(cls):
        return cls.__name__

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        raise NotImplementedError("This function has to be implemented in the derived class!")

    @classmethod
    def _CheckAvailabilityTransferOptions(cls, transfer_options):
        avail_options = cls._GetListAvailableTransferOptions()
        for option_name in transfer_options.GetStringArray():
            if not option_name in avail_options:
                err_msg  = 'transfer option "{}" not recognized for "{}"!\n'.format(option_name, cls._ClassName())
                err_msg += 'Available options: "{}"'.format('", "'.join(avail_options))
                raise Exception(err_msg)

    @classmethod
    def _GetDefaultSettings(cls):
        return KM.Parameters("""{
            "type"       : "UNSPECIFIED",
            "echo_level" : 0
        }""")
