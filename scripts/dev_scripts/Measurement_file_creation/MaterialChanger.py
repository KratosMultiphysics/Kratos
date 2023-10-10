from abc import ABC
from abc import abstractclassmethod

import KratosMultiphysics as Kratos


class MaterialChanger(ABC):

    """_summary_
    Interface for classes that are used in the process of generating artificial measurement data.
    """
    def __init__(self):
        pass

    @abstractclassmethod
    def adjust_material_of_model_part(self, model_part: Kratos.ModelPart) -> None:
        pass
