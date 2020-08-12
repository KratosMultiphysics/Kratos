# pyKratos imports
from .data_value_container import DataValueContainer

class Element(DataValueContainer):

    def __init__(self, elem_id, nodes):
        super().__init__()
        self.Id = elem_id
        self.__nodes = nodes
        self.__variables = {}


    def GetNode(self, node_index):
        return self.__nodes[node_index]

    def GetNodes(self):
        return self.__nodes

    def Initialize(self):
        pass

    def __str__(self):
        return  "Element #{} with {}".format(self.Id, self.__variables)

