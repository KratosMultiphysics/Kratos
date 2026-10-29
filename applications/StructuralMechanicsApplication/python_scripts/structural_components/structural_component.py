import numpy as np
from abc import ABC, abstractmethod

class StructuralComponent(ABC):

    def __init__(self, sub_model_part, metadata):
        self.sub_model_part = sub_model_part
        self.metadata = metadata
        self.analysis_results = {}

    @property
    def name(self) -> str:
        return self.sub_model_part.Name

    @classmethod
    @abstractmethod
    def FromKratosParametersObject(self, sub_model_part, data):
        """Instantiates the component from its entry in the configuration file.

        Args:
            sub_model_part (Kratos.SubModelPart): SubModelPart from a Kratos ModelPart
            data (Kratos.Parameters): Kratos.Parameters object with information from the configuration file
        """
        pass
    
    @abstractmethod
    def Initialize(self):
        pass

    @abstractmethod
    def PrepareAnalysis(self):
        """Extract data that methods need after FE solution.
        """
        pass
