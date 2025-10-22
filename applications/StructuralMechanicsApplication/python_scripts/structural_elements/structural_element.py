import numpy as np
from abc import ABC, abstractmethod

class StructuralElement(ABC):

    def __init__(self, 
                 sub_model_part, 
                 boundary_conditions: list[float], 
                 analysis_methods: list[str]):
        self.sub_model_part = sub_model_part
        self.boundary_conditions: list[float] = boundary_conditions
        self.analysis_methods: list[str] = analysis_methods
        #self.sub_model_part.AddNodes([self.corner_node_x_id, self.corner_node_y_id, self.origin_node_id])
        #self.InitializeNodeVectors()
        #self.SetLocalCoordinateSystem()
         

    @classmethod
    @abstractmethod
    def FromDict(self, sub_model_part, data):
        pass

    @classmethod
    @abstractmethod
    def FromKratosParametersObject(self, sub_model_part, data):
        """This classmethod instantiates a StructuralElement Object from a Kratos.Prameters Object 

        Args:
            sub_model_part (Kratos.SubModelPart): SubModelPart from a Kratos ModelPart
            data (Kratos.Parameters): Kratos.Parameters object with information from the configuration file
        """
        pass
    
    @abstractmethod
    def GetMaterialData(self):
        """Reads the material properties from the Kratos ModelPart and assigns them to the attributes of a structural element.
        """
        pass

    @abstractmethod
    def ComputeLoad(self):
        """Computes the loading of the structural element from the FE-Simulation results (e.g. from stresses).
        """
        pass

    @abstractmethod
    def RunHandbookMethods(self):
        """This method triggers the execution of the handbook methods that are defined in the configuration file.
        """
        pass

    @abstractmethod
    def ComputeMeasurements(self):
        """Compute the measurements of a structural element.
        """
        pass

    #def InitializeNodeVectors(self):
    #    """Initializes the vectors that point to the nodes defined in the configuration file. These vectors are used in other methods to compute
    #    attributes of a structural element (e.g. setting up the local coordinate system of a structural element)
    #    """
    #    self.origin_node_vector = np.array([self.sub_model_part.GetNode(self.origin_node_id).X, 
    #                                        self.sub_model_part.GetNode(self.origin_node_id).Y, 
    #                                        self.sub_model_part.GetNode(self.origin_node_id).Z])
    #    
    #    self.x_node_vector = np.array([self.sub_model_part.GetNode(self.corner_node_x_id).X,
    #                                   self.sub_model_part.GetNode(self.corner_node_x_id).Y,
    #                                   self.sub_model_part.GetNode(self.corner_node_x_id).Z])
    #    
    #    self.y_node_vector = np.array([self.sub_model_part.GetNode(self.corner_node_y_id).X,
    #                                   self.sub_model_part.GetNode(self.corner_node_y_id).Y,
    #                                   self.sub_model_part.GetNode(self.corner_node_y_id).Z])
        
    #def SetLocalCoordinateSystem(self):
    #    """This method sets up the local coordinate system of a structural element. It checks if the x- and y-axis are orthogonal. If that is not the case,
    #    a y-axis that is orthogonal to the x-axis is constructed. Afterwards a z-axis is constructued via cross-product.
    #    """
    #    tolerance = 0
    #    self.x_axis_base_vector = (self.x_node_vector - self.origin_node_vector)/np.linalg.norm(self.x_node_vector - self.origin_node_vector)
    #    self.y_axis_base_vector = (self.y_node_vector - self.origin_node_vector)/np.linalg.norm(self.y_node_vector - self.origin_node_vector)
    #    if abs(np.dot(self.x_axis_base_vector, self.y_axis_base_vector)) > tolerance:
    #        u1 = self.x_axis_base_vector
    #        u2 = (self.y_axis_base_vector - np.dot(self.y_axis_base_vector, u1)*u1)/np.linalg.norm(self.y_axis_base_vector - np.dot(self.y_axis_base_vector, u1)*u1)
    #        self.skew_angle = self.ComputeAngle(self.y_axis_base_vector, u2)
    #        self.y_axis_base_vector = u2
    #    self.z_axis_base_vector = np.cross(self.x_axis_base_vector, self.y_axis_base_vector)

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