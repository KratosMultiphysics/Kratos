from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

class KratosGlobals:

    def __init__(self, ThisKernel, ThisCaller, ApplicationsRoot):
        self.__dict__["Kernel"] = ThisKernel
        self.__dict__["RequestedApplications"] = dict()
        self.__dict__["AuthorizedCaller"] = ThisCaller
        self.__dict__["ApplicationsRoot"] = ApplicationsRoot
        self.__dict__["ApplicationsInterfaceIsDeprecated"] = False

    def __setattr__(self, name, value):
        if name in self.__dict__:
            # self.__dict__[name] = value
            print("Ignoring request to set KratosGlobals attribute", name)
        else:
            print("Ignoring request to set unknown KratosGlobals attribute:", name)

    def echo(self):
        print("Kernel:", self.Kernel)
        print("RequestedApplications:", self.RequestedApplications)
        print("Main Python script:", self.AuthorizedCaller)
        print("Kratos Applications base folder:", self.ApplicationsRoot)
        return
