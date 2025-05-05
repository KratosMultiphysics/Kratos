import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods import HandbookMethods as HM
import csv
from abc import ABC, abstractmethod

class StructuralElement(ABC):

    def __init__(self, sub_model_part, origin_node_id: int, corner_node_x_id: int, corner_node_y_id: int, analysis_methods: list[str]):
        self.sub_model_part = sub_model_part
        self.origin_node_id: int = origin_node_id
        self.corner_node_x_id: int = corner_node_x_id
        self.corner_node_y_id: int = corner_node_y_id
        self.analysis_methods: list[str] = analysis_methods
        self.sub_model_part.AddNodes([self.corner_node_x_id, self.corner_node_y_id, self.origin_node_id])
        self.InitializeNodeVectors()
        self.SetLocalCoordinateSystem()
         

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

    def InitializeNodeVectors(self):
        """Initializes the vectors that point to the nodes defined in the configuration file. These vectors are used in other methods to compute
        attributes of a structural element (e.g. setting up the local coordinate system of a structural element)
        """
        self.origin_node_vector = np.array([self.sub_model_part.GetNode(self.origin_node_id).X, 
                                            self.sub_model_part.GetNode(self.origin_node_id).Y, 
                                            self.sub_model_part.GetNode(self.origin_node_id).Z])
        
        self.x_node_vector = np.array([self.sub_model_part.GetNode(self.corner_node_x_id).X,
                                       self.sub_model_part.GetNode(self.corner_node_x_id).Y,
                                       self.sub_model_part.GetNode(self.corner_node_x_id).Z])
        
        self.y_node_vector = np.array([self.sub_model_part.GetNode(self.corner_node_y_id).X,
                                       self.sub_model_part.GetNode(self.corner_node_y_id).Y,
                                       self.sub_model_part.GetNode(self.corner_node_y_id).Z])
        
    def SetLocalCoordinateSystem(self):
        """This method sets up the local coordinate system of a structural element. It checks if the x- and y-axis are orthogonal. If that is not the case,
        a y-axis that is orthogonal to the x-axis is constructed. Afterwards a z-axis is constructued via cross-product.
        """
        tolerance = 0
        self.x_axis_base_vector = (self.x_node_vector - self.origin_node_vector)/np.linalg.norm(self.x_node_vector - self.origin_node_vector)
        self.y_axis_base_vector = (self.y_node_vector - self.origin_node_vector)/np.linalg.norm(self.y_node_vector - self.origin_node_vector)
        if abs(np.dot(self.x_axis_base_vector, self.y_axis_base_vector)) > tolerance:
            u1 = self.x_axis_base_vector
            u2 = (self.y_axis_base_vector - np.dot(self.y_axis_base_vector, u1)*u1)/np.linalg.norm(self.y_axis_base_vector - np.dot(self.y_axis_base_vector, u1)*u1)
            self.skew_angle = self.ComputeAngle(self.y_axis_base_vector, u2)
            self.y_axis_base_vector = u2
        self.z_axis_base_vector = np.cross(self.x_axis_base_vector, self.y_axis_base_vector)

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

class Panel(StructuralElement):

    def __init__(self, sub_model_part, origin_node_id: int, corner_node_x_id: int, corner_node_y_id: int, analysis_methods: list[str]):
        super().__init__(sub_model_part, origin_node_id, corner_node_x_id, corner_node_y_id, analysis_methods)
        self.GetMaterialData()
        self.ComputeLoad()
        self.ComputeMeasurements()
        self.RunHandbookMethods()

    @classmethod
    def FromDict(cls, sub_model_part, data):
        return cls(sub_model_part, 
            data["panel_origin_node"],
            data["corner_node_x"],
            data["corner_node_y"],
            data["analysis_methods"])
    
    @classmethod
    def FromKratosParametersObject(cls, sub_model_part, data):
        method_list = [data["analysis_methods"][i].GetString() for i in range(data["analysis_methods"].size())]
        return cls(sub_model_part, 
            data["panel_origin_node"].GetInt(),
            data["corner_node_x"].GetInt(),
            data["corner_node_y"].GetInt(),
            method_list)
    
    def ComputeMeasurements(self):
        #TODO: Attribute umbenennen
        #TODO: Length and width are also depending on the loadcase of the panel. This still has to be implemented here
        length_in_y = np.dot((self.y_node_vector - self.origin_node_vector), self.y_axis_base_vector)
        length_in_x = np.dot((self.x_node_vector - self.origin_node_vector), self.x_axis_base_vector)
        self.length = length_in_x
        self.width = length_in_y
        self.aspect_ratio = self.length/self.width
    
    def ComputeLoad(self):
        #TODO: Elementspannungen werden momentan in globalen Koordinaten ausgegeben. Sollten eigentlich in lokalen Koordinaten gegeben werden. 
        #TODO: Sobald Spannungen der Elemente in lokalen Koordinaten verwendet werden, muss die Transformation des Spannungstensors angepasst werden.
        #TODO: Elementspannungen werden an den Gaußpunkten berechnet. --> Spannungstensor aus der Mitte des Elements verwenden. Muss in Kratos implementiert werden?
        e1 = np.array([1, 0, 0])
        element_volumes = []
        element_xx_stresses_volume_product = []
        element_yy_stresses_volume_product = []
        for element in self.sub_model_part.Elements:
            element_stress = element.CalculateOnIntegrationPoints(SMA.SHELL_STRESS_MIDDLE_SURFACE, self.sub_model_part.ProcessInfo)
            phi = StructuralElement.ComputeAngle(self.x_axis_base_vector, e1)
            panel_stress = StructuralElement.rotation_around_z(element_stress[0], phi)
            element_volumes.append(element.GetGeometry().Area() * element.Properties.GetValue(KratosMultiphysics.THICKNESS))
            panel_volume = sum(element_volumes)
            element_xx_stresses_volume_product.append(panel_stress[0, 0] * element.GetGeometry().Area() * element.Properties.GetValue(KratosMultiphysics.THICKNESS))
            element_yy_stresses_volume_product.append(panel_stress[1, 1] * element.GetGeometry().Area() * element.Properties.GetValue(KratosMultiphysics.THICKNESS))
        self.xx_panel_stress = sum(element_xx_stresses_volume_product)/panel_volume
        self.yy_panel_stress = sum(element_yy_stresses_volume_product)/panel_volume

    
    def GetMaterialData(self):
        #TODO: Muss noch angepasst werden. Momentan nur valide, wenn das Material isotrop und homogen ist.
        i = 0
        for element in self.sub_model_part.Elements:
            if i == 0:
                properties = element.Properties
                #.Info() returns "ElasticIsotropic3D", but what I want is "LinearElasticPlaneStress2DLaw"
                #Since the derived class (LinearElasticPlaneStress2DLaw) doesn´t override the .Info() method of
                #the base class, it inhertis the method. So to get the actual constitutive law: (see line below)
                #TODO: Evtl. .Info() Methode für die ConstitutiveLaws anpassen?
                self.cl = properties.GetValue(KratosMultiphysics.CONSTITUTIVE_LAW).__class__.__name__
                self.thickness = properties.GetValue(KratosMultiphysics.THICKNESS)
                self.E = properties.GetValue(KratosMultiphysics.YOUNG_MODULUS)
                self.nu = properties.GetValue(KratosMultiphysics.POISSON_RATIO)
                i += 1
            else:
                break

    def CheckBoundaryConditions(self):
        """This function doesn´t do anything yet. I was just playing around a bit. The boundary conditions for a structural element should be given by
        the user in the configuration file. I still need to add this to the config file.
        """
        for node in self.sub_model_part.Nodes:
            for dof in ["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", "ROTATION_X", "ROTATION_Y", "ROTATION_Z"]:
                bc = node.IsFixed(getattr(KratosMultiphysics, dof))
                print(f"Is {dof} of {node.Id} fixed:", bc)

    def RunHandbookMethods(self):
        """This function has also only been implemented for testing purposes so far.
        """
        #TODO: Von-Mises ist nur als Test implementiert. Muss noch sauber geschrieben und evtl. angepasst werden.
        methods = self.analysis_methods
        for method in methods:
            match method:
                case "panel_buckling":
                    pass
                case "von_mises":
                    self.RF_von_mises = {}
                    for element in self.sub_model_part.Elements:
                        von_mises_stress = element.CalculateOnIntegrationPoints(SMA.VON_MISES_STRESS_MIDDLE_SURFACE, self.sub_model_part.ProcessInfo)
                        self.RF_von_mises[element.Id] = HM.von_mises(von_mises_stress[0], 100)
                    with open("von_mises_test.csv", "w", newline="") as f:
                        w = csv.DictWriter(f, self.RF_von_mises.keys())
                        w.writeheader()
                        w.writerow(self.RF_von_mises)