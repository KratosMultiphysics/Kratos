from __future__ import print_function, absolute_import, division 
class Node:

    def __init__(self, Id, coordinates):
        self.variables = []
        self.coordinates = coordinates
        self.Id = Id
        self.X = coordinates[0]
        self.Y = coordinates[1]
        self.Z = coordinates[2]

        self.X0 = coordinates[0]
        self.Y0 = coordinates[1]
        self.Z0 = coordinates[2]

    def SetBufferSize(self, buffer_size):
        for i in range(0, buffer_size):
            self.variables.append(dict())

    def AddVariable(self, variable_name):
        for i in range(0, len(self.variables)):
            self.variables[i][variable_name] = 0
        self.var_is_fixed[variable_name] = False

    def AdvanceInTime(self):
        for i in range(len(self.variables)-1,0,-1):
            for key in list(self.variables[i].keys()):
                self.variables[i][key] = self.variables[i - 1][key]

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

    def __str__(self):
        return  "Node #{0} with {1}".format(self.Id, self.variables)

