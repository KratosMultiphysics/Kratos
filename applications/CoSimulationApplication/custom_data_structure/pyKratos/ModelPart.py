from __future__ import print_function, absolute_import, division 
from .Node import *
from .Variables import *


class ModelPart:
    def __init__(self, name="default", buffer_size=1):
        self.NodesMap = {}
        self.Properties = {}
        self.ElementsMap = {}
        self.buffer_size = buffer_size
        self.solution_step_variables = []
        if "." in name:
            RuntimeError("Name of the ModelPart cannot contain a . (dot).")
        if name == "":
            RuntimeError("Name of ModelPart cannot be empty.")

        self.Name = name
        self.ProcessInfo = {TIME: 0.0, DELTA_TIME: 0.0}

        self.Nodes = list(self.NodesMap.values())

    def AddNodalSolutionStepVariable(self, variable):
        if isinstance(variable, list):  # For vector variables
            self.solution_step_variables.append(variable[0])
            self.solution_step_variables.append(variable[1])
            self.solution_step_variables.append(variable[2])
        else:
            self.solution_step_variables.append(variable)
        # To do: add variables to existing Nodes

    # Function to create a list of nodes based on a dictionary and give it to the model part
    def AddNodes(self, dict_of_nodes):
        for node_id, node_coords in list(dict_of_nodes.items()):
            if node_id in list(self.NodesMap.keys()):
                error_string = "Node with ID = " + str(node_id) + " already exists."
                raise Exception(error_string)
            else:
                node = Node(node_id, node_coords)
                node.SetBufferSize(self.buffer_size)
                for var in self.solution_step_variables:
                    node.AddVariable(var)
                self.NodesMap.update({node_id: node})
        self.Nodes = list(self.NodesMap.values())

    def AddNode(self, node):
        if type(node) == Node:
            self.NodesMap[node.Id] = node
            self.Nodes = list(self.NodesMap.values())
        else:
            RuntimeError("Adding a non-Node type object as a Node.")

    def CreateNewNode(self, node_id, x, y, z):
        if node_id in list(self.NodesMap.keys()):
            error_string = "Node with ID = " + str(node_id) + " already exists."
            raise Exception(error_string)
        else:
            node = Node(node_id, [x, y, z])
            node.SetBufferSize(self.buffer_size)
            for var in self.solution_step_variables:
                node.AddVariable(var)
            self.NodesMap.update({node_id: node})
        self.Nodes = list(self.NodesMap.values())

    def NumberOfNodes(self):
        return len(self.NodesMap)

    def NumberOfElements(self):
        return len(self.ElementsMap)

    def AddProperties(self, dict_of_properties):
        self.Properties.update(dict_of_properties)

    def AddElements(self, dict_of_elements):
        pass

    def GetElement(self):
        pass

    def GetNode(self, id):
        pass

    def RemoveElement(self):
        pass

    def RemoveNode(self):
        pass

    def HasNodalSolutionStepVariable(self):
        pass

    def CloneTimeStep(self, time):
        pass

    def WriteMesh(self):
        pass

    def Check(self):
        pass

    def __str__(self):
        return "ModelPart:\n    Number of Nodes: {0}\n    Nunber of Elements: {1}".format(len(self.NodesMap), len(self.ElementsMap))
