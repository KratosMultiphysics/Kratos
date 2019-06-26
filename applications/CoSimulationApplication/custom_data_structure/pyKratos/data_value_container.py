from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Other imports
from copy import deepcopy

class DataValueContainer(object):
    def __init__(self):
        self.__non_hist_variables = {}

    def SetValue(self, variable, value):
        # overwrite existing value or add new one
        # copy the value to match the behavior of Kratos (when used in python)
        # and to avoid unwanted references to wrong objects (for non-scalar types)
        self.__non_hist_variables[variable] = deepcopy(value)

    def GetValue(self, variable):
        if not variable in self.__non_hist_variables:
            # allocate this variable if it does not yet exist
            # this matches the Kratos behavior
            self.__non_hist_variables[variable] = variable.Zero()

        return self.__non_hist_variables[variable]

    def Has(self, variable):
        return variable in self.__non_hist_variables