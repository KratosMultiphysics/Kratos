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
        """Run handbook stress methods.
        """
        pass

    @staticmethod
    def ComputeAngle(a: np.ndarray, b: np.ndarray) -> float:
        """Computes the angle between two vectors a and b

        Args:
            a (np.ndarray) 
            b (np.ndarray)

        Returns:
            float: Angle between vectors a and b
        """
        v1 = np.array(a)
        v2 = np.array(b)
        v1_norm = np.linalg.norm(v1)
        v2_norm = np.linalg.norm(v2)
        phi = np.arccos(np.dot(v1, v2)/(v1_norm*v2_norm))
        return phi
    
    @staticmethod
    def rotation_around_z(sigma: np.ndarray, angle: float) -> np.ndarray:
        """Method I use to transform the stress tensor to the local structural element coordinate system.
        Right now this works for the use case but I think this still needs to be adjusted for future structural elements...

        Args:
            sigma (np.ndarray): Stress tensor in global coordinate system
            angle (float)

        Returns:
            np.ndarray: Stress tensor in local panel coordinate system
        """
        #TODO
        Q = np.array([[np.cos(angle), np.sin(angle), 0], [-np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        Q_T = np.transpose(Q)
        sigma_new = np.dot(Q,np.dot(sigma, Q_T))
        return sigma_new
