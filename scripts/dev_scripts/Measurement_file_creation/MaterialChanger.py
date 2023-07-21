from abc import ABC
from abc import abstractclassmethod

import KratosMultiphysics as Kratos


class MaterialChanger(ABC):

    def __init__(self):
        pass

    @abstractclassmethod
    def adjust_material_of_model_part(self, model_part: Kratos.ModelPart) -> None:
        pass
