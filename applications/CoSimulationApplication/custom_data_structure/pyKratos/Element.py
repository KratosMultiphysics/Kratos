from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

class Element(object):

    def __init__(self, elem_id, nodes):
        self.Id = elem_id
        self._nodes = nodes
        self._variables = {}

    def GetNode(self, node_index):
        return self._nodes[node_index]

    def GetNodes(self):
        return self._nodes

    def SetValue(self, variable, value):
        # overwrite existing value or add new one
        self._variables[variable] = value

    def GetValue(self, variable):
        if not variable in self._variables:
            # allocate this variable if it does not yet exist
            # this matches the Kratos behavior
            self._variables[variable] = variable.Zero()

        return self._variables[variable]

    def Has(self, variable):
        return variable in self._variables

    def Initialize(self):
        pass

    def __str__(self):
        return  "Element #{} with {}".format(self.Id, self._variables)

