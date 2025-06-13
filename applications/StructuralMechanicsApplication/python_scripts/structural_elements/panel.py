import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods import StrengthMethods, StabilityMethods
from KratosMultiphysics.StructuralMechanicsApplication.structural_elements.structural_element import StructuralElement
import csv

class Panel(StructuralElement):

    def __init__(self, 
                 sub_model_part, 
                 origin_node_id: int, 
                 corner_node_x_id: int, 
                 corner_node_y_id: int,
                 boundary_conditions: list[float],
                  analysis_methods: list[str]):
        super().__init__(sub_model_part, origin_node_id, corner_node_x_id, corner_node_y_id, boundary_conditions, analysis_methods)
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
        """This methods instantiates a panel object from given kratos parameters.

        Args:
            sub_model_part (_type_): SubModelPart defined in the *.mdpa file
            data (_type_): Data from the configuration file (*.json)
        """
        method_list = [data["analysis_methods"][i].GetString() for i in range(data["analysis_methods"].size())]
        boundary_conditions = [data["boundary_conditions"][i].GetDouble() for i in range(data["boundary_conditions"].size())]
        return cls(sub_model_part, 
            data["panel_origin_node"].GetInt(),
            data["corner_node_x"].GetInt(),
            data["corner_node_y"].GetInt(),
            boundary_conditions,
            method_list)
    
    def ComputeMeasurements(self) -> None:
        """This function calculates the measurements of the panel in the user defined x- and y-direction.
        """
        length_in_y = np.dot((self.y_node_vector - self.origin_node_vector), self.y_axis_base_vector)
        length_in_x = np.dot((self.x_node_vector - self.origin_node_vector), self.x_axis_base_vector)
        self.x_measurement: float = length_in_x
        self.y_measurement: float = length_in_y
    
    def ComputeLoad(self) -> None:
        #TODO: Elementspannungen werden momentan in globalen Koordinaten ausgegeben. Sollten eigentlich in lokalen Koordinaten gegeben werden. 
        #TODO: Sobald Spannungen der Elemente in lokalen Koordinaten verwendet werden, muss die Transformation des Spannungstensors angepasst werden.
        #TODO: Elementspannungen werden an den Gaußpunkten berechnet. --> Spannungstensor aus der Mitte des Elements verwenden. Muss in Kratos implementiert werden?
        e1 = np.array([1, 0, 0])
        element_volumes = []
        element_xx_stresses_volume_product = []
        element_yy_stresses_volume_product = []
        for element in self.sub_model_part.Elements:
            element_stress = element.CalculateOnIntegrationPoints(SMA.SHELL_STRESS_MIDDLE_SURFACE, self.sub_model_part.ProcessInfo)
            element_stress_global = element.CalculateOnIntegrationPoints(SMA.SHELL_STRESS_MIDDLE_SURFACE_GLOBAL, self.sub_model_part.ProcessInfo)
            local_x_axis = element.CalculateOnIntegrationPoints(KratosMultiphysics.LOCAL_AXIS_1, self.sub_model_part.ProcessInfo)[0]
            local_y_axis = element.CalculateOnIntegrationPoints(KratosMultiphysics.LOCAL_AXIS_2, self.sub_model_part.ProcessInfo)[0]
            phi = StructuralElement.ComputeAngle(self.x_axis_base_vector, e1)
            panel_stress = StructuralElement.rotation_around_z(element_stress[0], phi)
            element_volumes.append(element.GetGeometry().Area() * element.Properties.GetValue(KratosMultiphysics.THICKNESS))
            panel_volume = sum(element_volumes)
            element_xx_stresses_volume_product.append(panel_stress[0, 0] * element.GetGeometry().Area() * element.Properties.GetValue(KratosMultiphysics.THICKNESS))
            element_yy_stresses_volume_product.append(panel_stress[1, 1] * element.GetGeometry().Area() * element.Properties.GetValue(KratosMultiphysics.THICKNESS))
            #TODO: Delete the if condition and the Logger statement. Was used for debugging purposes.
            if element.Id == 21 or element.Id == 2:
                KratosMultiphysics.Logger.PrintInfo(f"Angle: {phi}; \n Element {element.Id} Stress: {element_stress[0]}; \n Element 1-axis: {local_x_axis}; \n Element 2-axis: {local_y_axis}; \nrotated stress: {panel_stress}, \n element stress global: {element_stress_global[0]}")
        self.xx_panel_stress: float = sum(element_xx_stresses_volume_product)/panel_volume
        self.yy_panel_stress: float = sum(element_yy_stresses_volume_product)/panel_volume

    
    def GetMaterialData(self) -> None:
        #TODO: Muss noch angepasst werden. Momentan nur valide, wenn das Material isotrop und homogen ist.
        i = 0
        for element in self.sub_model_part.Elements:
            if i == 0:
                properties = element.Properties
                #.Info() returns "ElasticIsotropic3D", but what I want is "LinearElasticPlaneStress2DLaw"
                #Since the derived class (LinearElasticPlaneStress2DLaw) doesn´t override the .Info() method of
                #the base class, it inhertis the method. So to get the actual constitutive law: (see line below)
                #TODO: Evtl. .Info() Methode für die ConstitutiveLaws anpassen?
                self.cl: str = properties.GetValue(KratosMultiphysics.CONSTITUTIVE_LAW).__class__.__name__
                self.thickness: float = properties.GetValue(KratosMultiphysics.THICKNESS)
                self.E: float = properties.GetValue(KratosMultiphysics.YOUNG_MODULUS)
                self.nu: float = properties.GetValue(KratosMultiphysics.POISSON_RATIO)
                if properties.Has(SMA.SHELL_ORTHOTROPIC_LAYERS):
                    self.orthotropic_material = True
                    print("ORTHOTROPIC")
                else:
                    self.orthotropic_material = False
                    print("ISOTROPIC")
                
                i += 1
            else:
                break

    def CheckBoundaryConditions(self):
        """This function doesn´t do anything yet. I was just playing around a bit. The boundary conditions for a structural element should be given by
        the user in the configuration file. I still need to add this to the config file.
        """
        for node in self.sub_model_part.Nodes:
            for dof in ["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", "ROTATION_X", "ROTATION_Y", "ROTATION_Z"]:
                bc = node.IsFixed(KratosMultiphysics.KratosGlobals.GetVariable(dof))
                print(f"Is {dof} of {node.Id} fixed:", bc)

    def PrepareStabilityAnalysis(self) -> None:
        """This functions checks if the panel is in uniaxial, biaxial or no compression and assigns the panel measurements and stresses accordingly.
        Conventionally compression stresses are assumed to be positive for buckling methods.
        """
        if self.xx_panel_stress < 0. and self.yy_panel_stress > 0.:
            self.buckling_mode: str = 'uniaxial'
            self.a: float = self.x_measurement
            self.b: float = self.y_measurement
            self.uniaxial_buckling_stress: float = -1*self.xx_panel_stress

        elif self.xx_panel_stress > 0. and self.yy_panel_stress < 0.:
            self.buckling_mode = 'uniaxial'
            self.a = self.y_measurement
            self.b = self.x_measurement
            self.uniaxial_buckling_stress = -1*self.yy_panel_stress

        elif self.xx_panel_stress < 0. and self.yy_panel_stress < 0.:
            self.buckling_mode = 'biaxial'
            self.beta: float = self.yy_panel_stress/self.xx_panel_stress
            if self.beta > 1.:
                self.beta = 1/self.beta
                self.a = self.y_measurement
                self.b = self.x_measurement
                self.biaxial_buckling_stress: list[float] = [-1*self.yy_panel_stress, -1*self.xx_panel_stress]
            else:
                self.a = self.x_measurement
                self.b = self.y_measurement
                self.biaxial_buckling_stress = [-1*self.xx_panel_stress, -1*self.yy_panel_stress]
                #KratosMultiphysics.Logger.PrintInfo(self.sub_model_part.Name, "TEST KRATOS LOGGER")

        else:
            self.buckling_mode = "None"
            KratosMultiphysics.Logger.PrintInfo(f"{self.sub_model_part.Name} is in tension")
            #print(f"{self.sub_model_part.Name} is in tension (or no valid stress results available). No buckling calculation is needed.")


        self.aspect_ratio = self.a/self.b

    def RunHandbookMethods(self) -> None:
        """This function has also only been implemented for testing purposes so far.
        """
        #TODO: Von-Mises ist nur als Test implementiert. Muss noch sauber geschrieben und evtl. angepasst werden.
        methods = self.analysis_methods
        for method in methods:
            if method == "panel_buckling":
                self.PrepareStabilityAnalysis()
                if self.buckling_mode == 'uniaxial':
                    RF = StabilityMethods.UniaxialBuckling(self.E, self.nu, self.thickness, self.a, abs(self.uniaxial_buckling_stress))
                elif self.buckling_mode == 'biaxial':
                    #TODO: Methode für biaxial anpassen
                    RF = StabilityMethods.BiaxialBuckling(self.E, self.nu, self.a, self.b, self.thickness, self.beta, self.biaxial_buckling_stress)
                elif self.buckling_mode == None:
                    RF = 8888
                with open("buckling.txt", "a") as f:
                    f.write(f"{self.sub_model_part.Name} : \n \t Buckling Mode: {self.buckling_mode}; \n \t RF: {RF} \n \n")
                    f.close()

            elif method == "von_mises":
                self.RF_von_mises = {}
                for element in self.sub_model_part.Elements:
                    von_mises_stress = element.CalculateOnIntegrationPoints(SMA.VON_MISES_STRESS_MIDDLE_SURFACE, self.sub_model_part.ProcessInfo)
                    self.RF_von_mises[element.Id] = StrengthMethods.VonMises(235, von_mises_stress[0])
                with open("von_mises_test.csv", "w", newline="") as f:
                    w = csv.DictWriter(f, self.RF_von_mises.keys())
                    w.writeheader()
                    w.writerow(self.RF_von_mises)

            else:
                print(f"{method} is not a valid input for methods.")
                continue