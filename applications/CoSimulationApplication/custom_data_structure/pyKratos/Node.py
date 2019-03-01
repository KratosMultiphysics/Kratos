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
        self.var_is_fixed = {}

    def SetBufferSize(self, buffer_size):
        for i in range(0, buffer_size):
            self.variables.append(dict())

    def AddVariable(self, variable):
        for i in range(0, len(self.variables)):
            self.variables[i][variable] = 0
        self.var_is_fixed[variable] = False

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

    def GetSolutionStepValue(self, variable, step):
        if(isinstance(variable, list)):
            return [self.variables[step][variable[1]],  self.variables[step][variable[2]], self.variables[step][variable[3]]]
        else:
            return self.variables[step][variable]

    def SetSolutionStepValue(self, variable, step, value):
        if(isinstance(variable, list)):
            for i in range(1, len(variable)):
                if variable[i] in list(self.variables[step].keys()):
                        self.variables[step][variable[i]] = value[i-1]
                else:
                    raise Exception(
                        "trying to set an non-existing variable with name ",
                        variable,
                        " on node ",
                        self.Id)

    def __str__(self):
        return  "Node #{0} with {1}".format(self.Id, self.variables)

