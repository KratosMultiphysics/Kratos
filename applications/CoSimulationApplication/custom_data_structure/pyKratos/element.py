from __future__ import print_function, absolute_import, division 
class Element:

    def __init__(self, Id, nodes_vector):
        self.variables = []
        self.Id = Id
        self.nodes = nodes_vector

    def __getitem__(self, key):
        return self.nodes[key]

    def GetNumberOfNodes(self):
        return len(self.nodes)

    def AddVariable(self, variable_name):
        for i in range(0, len(self.variables)):
            self.variables[i][variable_name] = 0

    def SetValue(self):
        pass # For non historical values

    def GetValue(self):
        pass

    def Has(self):
        pass

    def HasSolutionStepValue(self):
        pass

    def GetSolutionStepValue(self, variable_name, step):
        return self.variables[step][variable_name]

    def SetSolutionStepValue(self, variable_name, step, value):
        if variable_name in list(self.variables[step].keys()):
            self.variables[step][variable_name] = value
        else:
            raise Exception(
                "trying to set an non-existing variable with name ",
                variable_name,
                " on node ",
                self.Id)

    def GetArea(self):
        pass

    def GetNodes(self):
        pass

    def GetNormal(self):
        pass

    def Initialize(self):
        pass

    def __str__(self):
        return  "Element #{0} with {1}".format(self.Id, self.variables)

