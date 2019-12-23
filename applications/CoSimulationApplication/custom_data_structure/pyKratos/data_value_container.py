from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# pyKratos imports
from .Variable import VariableComponent

# Other imports
from copy import deepcopy

class DataValueContainer(object):
    def __init__(self):
        self.__non_hist_variables = {}

    def SetValue(self, variable, value):
        self.__InitializeVariable(variable)

        if isinstance(variable, VariableComponent):
            # here no explicit copy is necessary bcs it is a double
            self.__non_hist_variables[variable.GetSourceVariable()][variable.GetComponentIndex()] = value
        else:
            # copy the value to match the behavior of Kratos (when used in python)
            # and to avoid unwanted references to wrong objects (for non-scalar types)
            self.__non_hist_variables[variable] = deepcopy(value)

    def GetValue(self, variable):
        self.__InitializeVariable(variable)

        if isinstance(variable, VariableComponent):
            return self.__non_hist_variables[variable.GetSourceVariable()][variable.GetComponentIndex()]
        else:
            return self.__non_hist_variables[variable]

    def Has(self, variable):
        if isinstance(variable, VariableComponent):
            return variable.GetSourceVariable() in self.__non_hist_variables
        else:
            return variable in self.__non_hist_variables

    def __getitem__(self, variable):
        return self.GetValue(variable)

    def __setitem__(self, variable, value):
        return self.SetValue(variable, value)

    def __InitializeVariable(self, variable):
        # allocate this variable if it does not yet exist
        # this matches the Kratos behavior
        if isinstance(variable, VariableComponent):
            if not variable.GetSourceVariable() in self.__non_hist_variables:
                # allocate all components and the source-variable at the same time
                self.__non_hist_variables[variable.GetSourceVariable()] = variable.GetSourceVariable().Zero()
        else:
            if not variable in self.__non_hist_variables:
                self.__non_hist_variables[variable] = variable.Zero()
