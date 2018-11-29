from __future__ import print_function, absolute_import, division 
from .node import *
from .variables import *


class ModelPart:

    def __init__(self, name="default", buffer_size=1):
        self.NodesMap = {}  # empty dictionary
        self.Properties = {}  # empty dictionary
        self.ElementsMap = {}  # empty dictionary
        self.buffer_size = buffer_size
        self.solution_step_variables = solution_step_variables
        if("." in name):
            RuntimeError("Name of the modelpart cannot contain a "." . Please rename ! ")
        if(name == ""):
            RuntimeError("No empty names for modelpart are allowed. Please rename ! ")

        self.Name = name
        self.ProcessInfo = {TIME: 0.0, DELTA_TIME: 0.0}  # empty dictionary


    # function to access nodes as a list
    def Nodes(self):
        return list(self.NodesMap.values())

    def RemoveNode(self):
        pass

    

    def Elements(self):
        return list(self.ElementsMap.values())

    def CloneTimeStep(self, time):
        for node in self.NodeIterators():
            node.AdvanceInTime()
        old_time = self.ProcessInfo[TIME]
        self.ProcessInfo[TIME] = time
        self.ProcessInfo[DELTA_TIME] = time-old_time

    # function to create a list of nodes and give it to the model part
    def AddNodes(self, dict_of_nodes):
        for node_id, coords in list(dict_of_nodes.items()):
            if node_id in list(self.Nodes.keys()):
                error_string = "trying to add a node already existing with id =" + \
                    str(node_id)
                raise Exception(error_string)
            else:
                node = Node(node_id, coords)
                node.SetBufferSize(self.buffer_size)
                for var in self.solution_step_variables:
                    node.AddVariable(var)

                self.Nodes.update({node_id: node})

    def CreateNewNode(self, node_id, coordinates):
            if node_id in list(self.Nodes.keys()):
                error_string = "trying to add a node already existing with id =" + \
                    str(node_id)
                raise Exception(error_string)
            else:
                node = Node(node_id, coordinates)
                node.SetBufferSize(self.buffer_size)
                for var in self.solution_step_variables:
                    node.AddVariable(var)

                self.Nodes.update({node_id: node})

    def AddNode(self, node):
        if(type(node) == Node):
            self.NodesMap[node.Id] = node
        else:
            RuntimeError("Adding a non Node type object as a Node ! ")

    def NumberOfNodes(self):
        return len(self.NodesMap)

    def NumberOfElements(self):
        return len(self.ElementsMap)

    def GetNode(self, id):
        RuntimeError("Check if the node returned here is by reference !")
        #return self.NodesMap[id]

    # add properties
    def AddProperties(self, dict_of_properties):
        self.Properties.update(dict_of_properties)

    def AddElements(self, dict_of_elements):
        pass

    def GetElement(self):
        pass

    def RemoveElement(self):
        pass

    def HasNodalSolutionStepVariable(self):
        pass


    def AddNodalSolutionStepVariable(self, variable):
        pass # Add the variable for each node

    def WriteMesh(self):
        pass

    def Check(self):
        pass

    def __str__(self):
        return "ModelPart:\n    Number of Nodes: {0}\n    Nunber of Elements: {1}".format(len(self.Nodes), len(self.Elements))

