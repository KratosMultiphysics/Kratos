# pyKratos imports
from .ModelPart import ModelPart

class Model(dict):
    def GetModelPart(self, mp_name):
        # TODO this should get the ModelParts through splitting the names, same as in Kratos
        return self.__dict__[mp_name]

    def __getitem__(self, mp_name):
        return self.GetModelPart(mp_name)

    def CreateModelPart(self, mp_name="default", buffer_size=1):
        if mp_name in self.__dict__.keys():
            NameError("The modelpart with name ", mp_name, " already exist !")

        self.__dict__[mp_name] = ModelPart(mp_name, buffer_size, True)
        return self.__dict__[mp_name]

    def HasModelPart(self, mp_name):
        return mp_name in self.__dict__.keys()

    def DeleteModelPart(self, mp_name):
        self.__dict__.pop(mp_name)

    def __str__(self):
        self_str = "Model which contains the following ModelParts:"
        for mp_name in self.__dict__.keys():
            self_str += '\n\t{}'.format(mp_name)
        return self_str