import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import numpy as np
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods import StrengthMethods, StabilityMethods
from KratosMultiphysics.StructuralMechanicsApplication.structural_elements.structural_element import StructuralElement
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

class Panel(StructuralElement):

    def __init__(self, 
                 sub_model_part, 
                 boundary_conditions: list[float],
                  analysis_methods: list[str]):
        super().__init__(sub_model_part, boundary_conditions, analysis_methods)
        self.DefinePanelGeometry()
        self.GetMaterialData()
        self.ComputeLoad()
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
        #try: 
        #    panel_origin_node = data["panel_origin_node"].GetInt()
        #    corner_node_x = data["corner_node_x"].GetInt()
        #    corner_node_y = data["corner_node_y"].GetInt()
        return cls(sub_model_part, 
            boundary_conditions,
            method_list)
    
    def DefinePanelGeometry(self) -> None:

        element_areas = []
        element_centroids = []
        for element in self.sub_model_part.Elements:
            element_areas.append(element.GetGeometry().Area())
            coords = [np.array([node.X, node.Y, node.Z]) for node in element.GetGeometry()]
            el_centroid = sum(coords)/len(coords)
            element_centroids.append((element, el_centroid))
        
        centroid_coords = np.array([centroid for _, centroid in element_centroids])
        panel_center = np.mean(centroid_coords, axis=0)
        min_distance = float("inf")
        center_elem = None
        for elem, centroid in element_centroids:
            print("CENTROID:", centroid)
            dist = np.linalg.norm(centroid - panel_center)
            if dist < min_distance:
                min_distance = dist
                center_elem = elem
        self.center_element = elem
        self.avg_element_area = np.mean(element_areas)
        self.total_elemental_area = sum(element_areas)
        points = self.DefinePointCloud()
        #Test mit selbst definiertes Parabel
        #xs = np.linspace(0, 20, 21)
        #y1 = np.ones(shape=(21, ))*5
        #y2 = np.ones(shape=(21, ))*7
        #y3 = np.ones(shape=(21, ))*9
        #y4 = np.ones(shape=(21, ))*11
        #y5 = np.ones(shape=(21, ))*3
        #y6 = np.ones(shape=(21, ))*1
        #y7 = np.ones(shape=(21, ))*(-1)
        #z1 = np.ones(shape=(21, ))
        #z2 = np.ones(shape=(21, ))*2
        #z3 = np.ones(shape=(21, ))*4
        #z4 = np.ones(shape=(21, ))*7
        #p1s = np.vstack((xs, y1, z1))
        #p2s = np.vstack((xs, y2, z2))
        #p3s = np.vstack((xs, y3, z3))
        #p4s = np.vstack((xs, y4, z4))
        #p5s = np.vstack((xs, y5, z2))
        #p6s = np.vstack((xs, y6, z3))
        #p7s = np.vstack((xs, y7, z4))
        #points = np.hstack((p1s, p2s, p3s, p4s, p5s, p6s, p7s)).T
        x, y, z, proj_points = self.BestFittedPlane(points)
        len_x, len_y = self.ComputeMeasurements(proj_points)
        self.x_axis_base_vector = x
        self.y_axis_base_vector = y
        self.z_axis_base_vector = z
        self.x_measurement = round(len_x,2)
        self.y_measurement = round(len_y,2)

    def DefinePointCloud(self) -> np.ndarray:
        """Defined a point cloud from the nodes that are part of the submodelpart.

        Returns:
            points.T (np.ndarray): Point cloud as an nx3 matrix, where n is the number of points
        """
        N = self.sub_model_part.NumberOfNodes()
        points = np.empty(shape=(3, N))
        
        for i, node in enumerate(self.sub_model_part.Nodes):
            points[0][i] = node.X
            points[1][i] = node.Y
            points[2][i] = node.Z

        return points.T


    def BestFittedPlane(self, points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Defines a unit vector basis for the panel using principal component analysis.
            x points into the direction of greates variance
            y points into the directon of middle variance
            z points into the direction of smallest variance (i.e. the normal vector of the plan)
            proj_points is the point cloud projected onto the xy-plane of the new vector basis

        Args:
            points (np.ndarray): Point cloud as nx3 matrix

        Returns:
            tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]: Tuple consists of the three basis vectors for the panel and the points that were projected onto the panel plane
        """
        #TODO: Check skewness of panel and give warning if skewness is to high
        print("NODE_CLOUD SHAPE: ", points.shape)
        centroid = np.mean(points, axis=0)
        centered_points = points - centroid
        U, S, Vt = np.linalg.svd(centered_points)

        # Test for oblique panel testmodel -------------------------------------------------------------------------
        #corner_nodes = [1, 15, 36, 23] #oblique panel
        #corner_nodes = [1, 5, 6, 4] #shell angle test
        #n1_vector = np.array([self.sub_model_part.GetNode(corner_nodes[0]).X,
        #                        self.sub_model_part.GetNode(corner_nodes[0]).Y,
        #                        self.sub_model_part.GetNode(corner_nodes[0]).Z])
        #n2_vector = np.array([self.sub_model_part.GetNode(corner_nodes[1]).X,
        #                        self.sub_model_part.GetNode(corner_nodes[1]).Y,
        #                        self.sub_model_part.GetNode(corner_nodes[1]).Z])
        #n3_vector = np.array([self.sub_model_part.GetNode(corner_nodes[2]).X,
        #                        self.sub_model_part.GetNode(corner_nodes[2]).Y,
        #                        self.sub_model_part.GetNode(corner_nodes[2]).Z])
        #n4_vector = np.array([self.sub_model_part.GetNode(corner_nodes[3]).X,
        #                        self.sub_model_part.GetNode(corner_nodes[3]).Y,
        #                        self.sub_model_part.GetNode(corner_nodes[3]).Z])
#
        #diag_vec1 = n3_vector - n1_vector
        #diag_vec2 = n4_vector - n2_vector
        #n1n2_vector = n2_vector - n1_vector
        #n3n4_vector = n4_vector - n3_vector
        #print("DIAG",diag_vec1)
        #print("N1N2", n1n2_vector)
#
        #beta = np.arccos(np.dot(diag_vec1, n1n2_vector)/(np.linalg.norm(diag_vec1)*np.linalg.norm(n1n2_vector)))
        #gamma = np.arccos(np.dot((-1)*diag_vec2, (-1)*n1n2_vector)/(np.linalg.norm(diag_vec2)*np.linalg.norm(n1n2_vector)))
        #alpha = (beta+gamma)/2
        #print(f"Beta: {np.degrees(beta)}; Gamma: {np.degrees(gamma)}, Alpha: {np.degrees(alpha)}")
        #R = np.array([[np.cos(alpha), -np.sin(alpha), 0],
        #    [np.sin(alpha),  np.cos(alpha), 0],
        #    [0,0,0]])
        #x = diag_vec1 @ R.T
        #x /= np.linalg.norm(x)
        #z = Vt.T[-1]
        #y = np.cross(x, z)
        #y /= np.linalg.norm(y)
#
        #principal_axes = np.array([x,y,z])
        # --------------------------------------------------------------------------------------------------------------------
        principal_axes = Vt.T
        proj_points = centered_points @ principal_axes
        x = principal_axes[0]
        y = principal_axes[1]
        z = principal_axes[-1]
        #check alignment with local coordinate system of center element
        cos_theta = np.dot(x, self.center_element.CalculateOnIntegrationPoints(KratosMultiphysics.LOCAL_AXIS_1, self.sub_model_part.ProcessInfo)[0])
        phi = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
        print(self.center_element.CalculateOnIntegrationPoints(KratosMultiphysics.LOCAL_AXIS_1, self.sub_model_part.ProcessInfo)[0])
        if ((abs(phi) < 70) or (abs(phi) > 90)) & ((abs(phi) < 160) or abs(phi) > 200):
            angle_rad = np.radians(phi)
            R = np.array([[np.cos(angle_rad), -np.sin(angle_rad), 0],
            [np.sin(angle_rad),  np.cos(angle_rad), 0],
            [0, 0, 1]])
            x = x @ R
            #x = np.append(x, 0)
            y = np.cross(x, z)
            principal_axes = np.array([x, y, z])
            print("PRINC1:", principal_axes)
            proj_points = centered_points @ principal_axes

        # --------------------------------------------------------------------------------------------------------------------------------
        #principal_axes_niklas = np.array([  [-0.13541024, 0.99078962, 0.],
        #                                    [-0.99078962, -0.13541024,  0.],
        #                                    [-8.107655, 12.8377, 0.]])
        #x = principal_axes_niklas[0]
        #y = principal_axes_niklas[1]
        #z = principal_axes_niklas[-1]
        #proj_points = centered_points @ principal_axes_niklas
        # --------------------------------------------------------------------------------------------------------------------------------

        return x, y, z, proj_points
    
    def ComputeMeasurements(self, proj_points: np.ndarray) -> tuple[float, float]:
        """Computes the measurements (length and width) of the panel by splitting up the projected point cloud into "band" with a certain width. The min and max values of x and y position
        are subtracted from eachother to get the length and width of the band, respectively. The slicing is done in each direction separately. The band lengths and widths are then
        averaged to get the average measurements of the panel. (Note: The thickness information is part of the element property)

        Args:
            proj_points (np.ndarray): Nodes of submodelpart projected onto panel plane

        Returns:
            tuple[float, float]: [measurement in x direction, measurement in y direction]
        """

        xs = proj_points[:, 0]
        ys = proj_points[:, 1]
        num_slices_x = max(10, min(100, len(xs) // 5))
        num_slices_y = max(10, min(100, len(ys) // 5))#int(len(ys)/2)

        min_x, max_x = xs.min(), xs.max()
        min_y, max_y = ys.min(), ys.max()

        #TODO: Überlegen, wie man am besten die num_slices und die tolerance definiert (vermutlich am Besten über die Elementgeometrie die Toleranz)
        slices_x = np.linspace(min_x, max_x, 5)
        slices_y = np.linspace(min_y, max_y, 5)
        xs_dist = np.sort(xs)[-1] - np.sort(xs)[0]
        ys_dist = np.sort(ys)[-1] - np.sort(ys)[0]
        tolerance_x = (max_x - min_x) / 3 #np.sqrt(self.avg_element_area)
        tolerance_y = (max_y - min_y) / 3 #np.sqrt(self.avg_element_area)
        y_distances = []
        x_distances = []
        print(slices_x)
        print("TOLERANCE X:", tolerance_x)
        print("TOLERANCE Y:", tolerance_y)
        #get width (y-measurement of panel)
        i = 0
        for x_slice in slices_x:
            mask = (xs <= x_slice + tolerance_x) & (xs >= x_slice- tolerance_x)
            if np.count_nonzero(mask) < 2:
                continue
            ys_in_slice = ys[mask]
            try:
                min_y, max_y = np.min(ys_in_slice), np.max(ys_in_slice)
            except ValueError:
                #print("No ys in slice")
                continue
            distance = max_y - min_y
            if distance != 0:
                y_distances.append(distance)
            if i == 0:
                i += 1
                plt.figure(figsize=(6,6))
                points_in_slice = proj_points[mask]
                plt.scatter(proj_points[:, 0], proj_points[:, 1], color='lightgray', label='All Points')
                plt.scatter(points_in_slice[:, 0], points_in_slice[:, 1], color='red', label='Slice Points')

                # Optional: connect min/max y points with a line
                min_y_idx = np.argmin(ys_in_slice)
                max_y_idx = np.argmax(ys_in_slice)
                min_x = xs[min_y_idx]
                max_x = xs[max_y_idx]
                p1 = points_in_slice[min_y_idx]
                p2 = points_in_slice[max_y_idx]
                plt.plot([min_x, max_x], [min_y, max_y], 'blue', linewidth=2, label='Slice Width')

                plt.title(f"Slice at x = {x_slice:.2f}")
                plt.xlabel("x (PCA-aligned)")
                plt.ylabel("y (PCA-aligned)")
                plt.legend()
                plt.axis('equal')
                plt.grid(True)
                plt.show()

                fig, axs = plt.subplots(1, 2, figsize=(12, 6))
                points = self.DefinePointCloud()
                # 1️⃣ Original coordinate system
                axs[0].scatter(points[:, 0], points[:, 1], color='gray', label='Original Points')
                axs[0].set_title('Original Coordinate System')
                axs[0].set_xlabel('X')
                axs[0].set_ylabel('Y')
                axs[0].axis('equal')
                axs[0].grid(True)
                axs[0].legend()

                # 2️⃣ PCA-aligned coordinate system
                axs[1].scatter(proj_points[:, 0], proj_points[:, 1], color='red', label='Projected Points')
                axs[1].set_title('PCA-Aligned Coordinate System')
                axs[1].set_xlabel('Principal X')
                axs[1].set_ylabel('Principal Y')
                axs[1].axis('equal')
                axs[1].grid(True)
                axs[1].legend()

                plt.tight_layout()
                plt.show()

        for y_slice in slices_y:
            mask = (ys <= y_slice + tolerance_y) & (ys >= y_slice- tolerance_y)
            xs_in_slice = xs[mask]
            try:
                min_x, max_x = np.min(xs_in_slice), np.max(xs_in_slice)
            except ValueError:
                continue
            distance = max_x - min_x
            if distance != 0:
                x_distances.append(distance)

        length_in_x = np.mean(x_distances)
        length_in_y = np.mean(y_distances)
        return length_in_x, length_in_y
    
    def ComputeLoad(self) -> None:
        """Computes the average stress in the panel x- and y-direction.
        """
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
            if element.Id == 21 or element.Id == 2 or element.Id == 406:
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
            self.a = self.x_measurement
            self.b = self.y_measurement
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
                elif self.buckling_mode == "None":
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