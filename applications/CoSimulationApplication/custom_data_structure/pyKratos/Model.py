from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# pyKratos imports
from .ModelPart import ModelPart

class Model(dict):
    def GetModelPart(self, name):
        # TODO this should get the ModelParts through splitting the names, same as in Kratos
        return self.__dict__[name]

    def CreateModelPart(self, name="default", buffer_size=1):
        if(name in self.__dict__.keys() ):
            RuntimeError("The modelpart with name ", name, " already exist !")

        self.__dict__[name] = ModelPart(name, buffer_size, True)
        return self.__dict__[name]

    def HasModelPart(self, name):
        return (name in self.__dict__.keys())

    def DeleteModelPart(self, name):
        self.__dict__.pop(name)
