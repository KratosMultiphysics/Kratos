from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

class Node(object):

    def __init__(self, node_id, x, y, z, hist_variables, buffer_size):
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

        # non-historical variables
        self._non_hist_variables = {}

    ### Methods related to historical variables ###
    def CloneSolutionStep(self):
        for var_vals in self.__solution_steps_nodal_data.values():
            var_vals[1:self.__buffer_size] = var_vals[0:self.__buffer_size-1]

    def SetValue(self, variable, value):
        # overwrite existing value or add new one
        self._non_hist_variables[variable] = value

    def GetValue(self, variable):
        if not variable in self._non_hist_variables:
            # allocate this variable if it does not yet exist
            # this matches the Kratos behavior
            self._non_hist_variables[variable] = variable.Zero()

        return self._non_hist_variables[variable]

    def Has(self, variable):
        return variable in self._non_hist_variables

    def SolutionStepsDataHas(self, variable):
        return variable in self.__hist_variables

    def GetSolutionStepValue(self, variable, step):
        self.__CheckHistoricalVariable(variable)
        self.__CheckBufferSize(step)

        return self.__solution_steps_nodal_data[variable][step]

    def SetSolutionStepValue(self, variable, step, value):
        self.__CheckHistoricalVariable(variable)
        self.__CheckBufferSize(step)

        self.__solution_steps_nodal_data[variable][step] = value


    def __CheckHistoricalVariable(self, variable):
        if not variable in self.__hist_variables:
            raise Exception('Trying to access historical variable "{}" which does not exist!'.format(variable))

    def __CheckBufferSize(self, step):
        if step+1 > self.__buffer_size:
            raise Exception('Insufficient buffer size: requested: {}; available: {}'.format(step+1, self.__buffer_size))

    def __InitializeSolutionStepsNodalData(self):
        for var in self.__hist_variables:
            zero_val = var.Zero()
            self.__solution_steps_nodal_data[var] = [zero_val for i in range(self.__buffer_size)]

    def __str__(self):
        return  "Node #{0} with {1}".format(self.Id, self.__hist_variables)

