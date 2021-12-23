class KratosGlobalsImpl(object):

    def __init__(self, ThisKernel, ApplicationsRoot):
        self.__dict__["Kernel"] = ThisKernel
        self.__dict__["ApplicationsRoot"] = ApplicationsRoot

    def __setattr__(self, name, value):
        if name in self.__dict__:
            print("Ignoring request to set KratosGlobals attribute", name)
        else:
            print("Ignoring request to set unknown KratosGlobals attribute:", name)

    def echo(self):
        print("Kernel:", self.Kernel)
        print("Kratos Applications base folder:", self.ApplicationsRoot)

    def GetFlag(self, FlagName):
        """ This method returns the flag with the given name

        Keyword arguments:
        self -- It signifies an instance of a class.
        FlagName -- The name of the flag to return
        """
        return self.Kernel.GetFlag(FlagName)

    def GetVariable(self, VarName):
        """ This method returns the variable with the given name

        Keyword arguments:
        self -- It signifies an instance of a class.
        VarName -- The name of the variable to return
        """
        kernel = self.Kernel

        if kernel.HasDoubleVariable(VarName):
            return kernel.GetDoubleVariable(VarName)
        elif kernel.HasArrayVariable(VarName):
            return kernel.GetArrayVariable(VarName)
        elif kernel.HasArray4Variable(VarName):
            return kernel.GetArray4Variable(VarName)
        elif kernel.HasArray6Variable(VarName):
            return kernel.GetArray6Variable(VarName)
        elif kernel.HasArray9Variable(VarName):
            return kernel.GetArray9Variable(VarName)
        elif kernel.HasBoolVariable(VarName):
            return kernel.GetBoolVariable(VarName)
        elif kernel.HasIntVariable(VarName):
            return kernel.GetIntVariable(VarName)
        elif kernel.HasUnsignedIntVariable(VarName):
            return kernel.GetUnsignedIntVariable(VarName)
        elif kernel.HasVectorVariable(VarName):
            return kernel.GetVectorVariable(VarName)
        elif kernel.HasMatrixVariable(VarName):
            return kernel.GetMatrixVariable(VarName)
        elif kernel.HasStringVariable(VarName):
            return kernel.GetStringVariable(VarName)
        elif kernel.HasFlagsVariable(VarName):
            return kernel.GetFlagsVariable(VarName)
        elif kernel.HasVariableData(VarName):
            raise ValueError("\nKernel.GetVariable() ERROR: Variable {0} is defined but is of unsupported type\n".format(VarName))
        else:
            raise ValueError("\nKernel.GetVariable() ERROR: Variable {0} is unknown. Maybe you need to import the application where it is defined?\n".format(VarName))

    def HasVariable(self, VarName):
        """ This method checks if a variable exists

        Keyword arguments:
        self -- It signifies an instance of a class.
        VarName -- The name of the variable to check
        """
        kernel = self.Kernel

        if kernel.HasDoubleVariable(VarName):
            return True
        elif kernel.HasArrayVariable(VarName):
            return True
        elif kernel.HasArray4Variable(VarName):
            return True
        elif kernel.HasArray6Variable(VarName):
            return True
        elif kernel.HasArray9Variable(VarName):
            return True
        elif kernel.HasBoolVariable(VarName):
            return True
        elif kernel.HasIntVariable(VarName):
            return True
        elif kernel.HasUnsignedIntVariable(VarName):
            return True
        elif kernel.HasVectorVariable(VarName):
            return True
        elif kernel.HasMatrixVariable(VarName):
            return True
        elif kernel.HasStringVariable(VarName):
            return True
        elif kernel.HasFlagsVariable(VarName):
            return True
        elif kernel.HasVariableData(VarName):
            return True
        else:
            return False

    def GetVariableType(self, VarName):
        """ This method checks the type of variable

        Keyword arguments:
        self -- It signifies an instance of a class.
        VarName -- The name of the variable to check
        """
        kernel = self.Kernel

        if kernel.HasBoolVariable(VarName):
            return "Bool"
        elif kernel.HasIntVariable(VarName):
            return "Integer"
        elif kernel.HasUnsignedIntVariable(VarName):
            return "Unsigned Integer"
        elif kernel.HasDoubleVariable(VarName):
            return "Double"
        elif kernel.HasArrayVariable(VarName):
            return "Array"
        elif kernel.HasArray4Variable(VarName):
            return "Array4"
        elif kernel.HasArray6Variable(VarName):
            return "Array6"
        elif kernel.HasArray9Variable(VarName):
            return "Array9"
        elif kernel.HasVectorVariable(VarName):
            return "Vector"
        elif kernel.HasMatrixVariable(VarName):
            return "Matrix"
        elif kernel.HasStringVariable(VarName):
            return "String"
        elif kernel.HasFlagsVariable(VarName):
            return "Flag"
        else:
            return "NONE"

    def GetConstitutiveLaw(self, ConstitutiveLawName):
        """ This method returns the constitutive law with the given name

        Keyword arguments:
        self -- It signifies an instance of a class.
        ConstitutiveLawName -- The name of the constitutive law to return
        """
        return self.Kernel.GetConstitutiveLaw(ConstitutiveLawName)

    def HasConstitutiveLaw(self, ConstitutiveLawName):
        """ This method checks if a constitutive law exists

        Keyword arguments:
        self -- It signifies an instance of a class.
        ConstitutiveLawName -- The name of the constitutive law to check
        """
        return self.Kernel.HasConstitutiveLaw(ConstitutiveLawName)
