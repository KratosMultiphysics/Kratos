from __future__ import print_function, absolute_import, division 

from .ModelPart import ModelPart


class Model(dict):
    def GetModelPart(self, name):
        return self.__dict__[name]

    def CreateModelPart(self, name):
        if name in self.__dict__.keys():
            RuntimeError("The modelpart with name ", name, " already exists")

        self.__dict__[name] = ModelPart(name)
        return self.__dict__[name]

    def HasModelPart(self, name):
        return name in self.__dict__.keys()

    def DeleteModelPart(self, name):
        self.__dict__.pop(name)
