from __future__ import print_function, absolute_import, division 
from .Node import *
from .Variables import *


class ModelPart:

    def __init__(self, name="default", buffer_size=1):
        self.Nodes = {}  # empty dictionary
        self.Properties = {}  # empty dictionary
        self.Elements = {}  # empty dictionary
        self.buffer_size = buffer_size
        self.solution_step_variables = []
        if("." in name):
            RuntimeError("Name of the modelpart cannot contain a . (dot) Please rename ! ")
        if(name == ""):
            RuntimeError("No empty names for modelpart are allowed. Please rename ! ")

        self.Name = name
        self.ProcessInfo = {TIME: 0.0, DELTA_TIME: 0.0}  # empty dictionary


    def RemoveNode(self):
        pass

    def AddNodalSolutionStepVariable(self, variable):
        self.solution_step_variables.append(variable)

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

    def CreateNewNode(self, node_id, x, y, z):
            if node_id in list(self.Nodes.keys()):
                error_string = "trying to add a node already existing with id =" + \
                    str(node_id)
                raise Exception(error_string)
            else:
                node = Node(node_id, [x,y,z])
                node.SetBufferSize(self.buffer_size)
                for var in self.solution_step_variables:
                    node.AddVariable(var)

                self.Nodes.update({node_id: node})

    def AddNode(self, node):
        if(type(node) == Node):
            self.Nodes[node.Id] = node
        else:
            RuntimeError("Adding a non Node type object as a Node ! ")

    def NumberOfNodes(self):
        return len(self.Nodes)

    def NumberOfElements(self):
        return len(self.Elements)

    def GetNode(self, id):
        RuntimeError("Check if the node returned here is by reference !")
        #return self.Nodes[id]

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

