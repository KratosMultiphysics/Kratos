from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# pyKratos imports
from .data_value_container import DataValueContainer

class Element(object):

    def __init__(self, elem_id, nodes):
        self.Id = elem_id
        self.__nodes = nodes
        self.__variables = {}

        # non-historical variables
        self.__data_value_container = DataValueContainer()


    def GetNode(self, node_index):
        return self.__nodes[node_index]

    def GetNodes(self):
        return self.__nodes


    ### Methods related to non-historical variables ###
    def SetValue(self, variable, value):
        self.__data_value_container.SetValue(variable, value)

    def GetValue(self, variable):
        return self.__data_value_container.GetValue(variable)

    def Has(self, variable):
        return self.__data_value_container.Has(variable)

    def __getitem__(self, variable):
        return self.GetValue(variable)

    def __setitem__(self, variable, value):
        return self.SetValue(variable, value)


    def Initialize(self):
        pass

    def __str__(self):
        return  "Element #{} with {}".format(self.Id, self.__variables)

