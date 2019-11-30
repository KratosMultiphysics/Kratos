from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# pyKratos imports
from .data_value_container import DataValueContainer
from .Variable import VariableComponent

# Other imports
from copy import deepcopy

class Node(DataValueContainer):

    def __init__(self, node_id, x, y, z, hist_variables, buffer_size):
        super(Node, self).__init__()
        self.Id = node_id

        # current position
        self.X = x
        self.Y = y
        self.Z = z

        # initial position
        self.X0 = x
        self.Y0 = y
        self.Z0 = z

        # historical variables
        self.__hist_variables = hist_variables
        self.__buffer_size = buffer_size
        self.__solution_steps_nodal_data = {}
        self.__InitializeSolutionStepsNodalData()


    ### Methods related to historical variables ###
    def CloneSolutionStep(self):
        for var_vals in self.__solution_steps_nodal_data.values():
            var_vals[1:self.__buffer_size] = var_vals[0:self.__buffer_size-1]

    def SolutionStepsDataHas(self, variable):
        if isinstance(variable, VariableComponent):
            return variable.GetSourceVariable() in self.__hist_variables
        else:
            return variable in self.__hist_variables

    def GetSolutionStepValue(self, variable, step=0):
        self.__CheckHistoricalVariable(variable)
        self.__CheckBufferSize(step)

        if isinstance(variable, VariableComponent):
            return self.__solution_steps_nodal_data[variable.GetSourceVariable()][step][variable.GetComponentIndex()]
        else:
            return self.__solution_steps_nodal_data[variable][step]

    def SetSolutionStepValue(self, variable, step, value):
        self.__CheckHistoricalVariable(variable)
        self.__CheckBufferSize(step)

        if isinstance(variable, VariableComponent):
            self.__solution_steps_nodal_data[variable.GetSourceVariable()][step][variable.GetComponentIndex()] = value
        else:
            # copy the value to match the behavior of Kratos (when used in python)
            # and to avoid unwanted references to wrong objects (for non-scalar types)
            self.__solution_steps_nodal_data[variable][step] = deepcopy(value)


    def __CheckHistoricalVariable(self, variable):
        if isinstance(variable, VariableComponent):
            var_to_check = variable.GetSourceVariable()
        else:
            var_to_check = variable
        if not var_to_check in self.__hist_variables:
            raise Exception('Trying to access historical variable "{}" which does not exist!'.format(variable))

    def __CheckBufferSize(self, step):
        if step+1 > self.__buffer_size:
            raise Exception('Insufficient buffer size: requested: {}; available: {}'.format(step+1, self.__buffer_size))

    def __InitializeSolutionStepsNodalData(self):
        for var in self.__hist_variables:
            self.__solution_steps_nodal_data[var] = [var.Zero() for i in range(self.__buffer_size)]

    def __str__(self):
        return  "Node #{0} with {1}".format(self.Id, self.__hist_variables)

