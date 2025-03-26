# Importing Kratos
import KratosMultiphysics
from KratosMultiphysics.ConvectionDiffusionApplication import python_solvers_wrapper_convection_diffusion as solver_wrapper

# Importing the base class
from KratosMultiphysics.analysis_stage import AnalysisStage

# Other imports
import sys
import numpy as np
import scipy.sparse.linalg as spla
from shapely.geometry import LineString, Polygon
import matplotlib
matplotlib.use('Agg')  # Use Agg for rendering plots to files
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from colorama import Fore, Style, init
from scipy.spatial import KDTree
from itertools import product

# Interpolation libraries
from scipy.interpolate import Rbf

# Enable LaTeX rendering in Matplotlib
plt.rcParams.update({
    "text.usetex": True,             # Use LaTeX for text
    "font.family": "serif",          # Use serif font
    "pgf.texsystem": "pdflatex",     # Use pdflatex
    "pgf.preamble": r"\usepackage{amsmath, amssymb}",  # Load additional LaTeX packages
})

class ExtendedGradientMethodConvectionDiffusionAnalysis(AnalysisStage):
    """
    This class is the main-script of the ExtendedGradientConvectionDiffusion method put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # Making sure that older cases still work by properly initalizing the parameters
        self.solver_settings = project_parameters["solver_settings"]
        
        self.model = model

        if not self.solver_settings.Has("domain_size"):
            KratosMultiphysics.Logger.PrintInfo("ExtendedGradientMethodConvectionDiffusionAnalysis", "Using the old way to pass the domain_size, this will be removed!")
            self.solver_settings.AddEmptyValue("domain_size")
            self.solver_settings["domain_size"].SetInt(project_parameters["problem_data"]["domain_size"].GetInt())

        super(ExtendedGradientMethodConvectionDiffusionAnalysis, self).__init__(model, project_parameters)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        return solver_wrapper.CreateSolverByParameters(self.model, self.project_parameters["solver_settings"],self.project_parameters["problem_data"]["parallel_type"].GetString())

    def _GetSimulationName(self):
        return "::[Extended Gradient Method Convection-Diffusion Simulation]:: "
    
    def RunSolutionLoop(self):
        # Get the unknown variable for the specific physical problem which is being solved
        unknown_variable_name = self.solver_settings["convection_diffusion_variables"]["unknown_variable"].GetString()
        unknown_variable = KratosMultiphysics.KratosGlobals.GetVariable(unknown_variable_name)
        
        # Write the output for the initial conditions
        self.OutputSolutionStep()

        self.old_solution_field = np.zeros(self._GetSolver().GetComputingModelPart().NumberOfNodes())
        for element in self._GetSolver().GetComputingModelPart().Elements:
            for node in element.GetNodes():
                self.old_solution_field[node.Id - 1] = node.GetSolutionStepValue(unknown_variable) + 1e-2

        self.numerical_solution_with_iterations = []
        self.numerical_solution_no_iterations = []
        self.exact_solution = []
        self.time_list = []
        self.number_of_iterations = []

        while self.KeepAdvancingSolutionLoop():
            # Initialize the variables for the algorithm
            epsilon = 1e-3
            is_not_converged = True
            self.iteration_number = 0
            self.maximum_iterations = 20

            # Initialize the old solution field
            # if ():
            #     self.old_solution_field = np.zeros(self._GetSolver().GetComputingModelPart().NumberOfNodes())
            #     for element in self._GetSolver().GetComputingModelPart().Elements:
            #         for node in element.GetNodes():
            #             self.old_solution_field[node.Id - 1] = node.GetSolutionStepValue(unknown_variable) + 1e-2

            self.time = self._AdvanceTime()
            self.InitializeSolutionStep()

            while (is_not_converged == True and self.iteration_number < self.maximum_iterations):
                self._GetSolver().Predict()
                self.ApplyDirichletBoundaryConditions()
                is_converged = self._GetSolver().SolveSolutionStep()
                solution_error = self.ComputeErrorAndUpdateOldSolutionField(self.iteration_number)
                print(f"{Fore.BLUE}Iteration # {self.iteration_number}{Style.RESET_ALL}")
                print(f"{Fore.RED}Error between iterations: {solution_error}{Style.RESET_ALL}")

                if (solution_error < epsilon or self.iteration_number == 19):
                    is_not_converged = False
                    print( f"{Fore.GREEN}*** CONVERGENCE ACHIEVED ***{Style.RESET_ALL}")

                self.iteration_number += 1

                if self.iteration_number == 1:
                    self.numerical_solution_no_iterations.append(self.model_part.GetNode(1926).GetSolutionStepValue(unknown_variable))

            self.number_of_iterations.append(self.iteration_number)

            for element in self._GetSolver().GetComputingModelPart().Elements:
                for node in element.GetNodes():
                    self.old_solution_field[node.Id - 1] = node.GetSolutionStepValue(unknown_variable) 

            # Plot the solution at a specific 
            self.time_list.append(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
            self.numerical_solution_with_iterations.append(self.model_part.GetNode(1926).GetSolutionStepValue(unknown_variable))
            self.exact_solution.append(self.BC_value(self.model_part.GetNode(1926).X, self.model_part.GetNode(1926).Y))

            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

        # First plot
        plt.figure()  # Create a new figure
        plt.plot(self.time_list, self.numerical_solution_with_iterations, color='r', label=r'Solution with iterations')
        plt.plot(self.time_list, self.numerical_solution_no_iterations, color='g', label=r'Solution without iterations')
        plt.plot(self.time_list, self.exact_solution, linestyle='--', color='b', label=r'Exact Solution')

        plt.xlabel(r'Time [s]')
        plt.ylabel(r'Solution Value')
        plt.grid(True, which="both", ls="--")
        plt.legend()
        plt.savefig("numerical_solution_plot.pdf", dpi=300)
        plt.show()

        # Second plot
        plt.figure()  # Create a new figure
        plt.plot(self.time_list, self.number_of_iterations, color='r', label=r'Number of iterations')

        plt.xlabel(r'Time [s]')
        plt.ylabel(r'Number of iterations for convergence $\mathcal{N}$')
        plt.grid(True, which="both", ls="--")
        plt.legend()
        plt.savefig("iterations_vs_time.pdf", dpi=300)
        plt.show()

    def BC_value(self, x, y):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        # if time <= 0.01:
        #     return np.sin(np.pi * x)* np.cos(np.pi * y)*(time/0.01)
        # else:
        return np.sin(np.pi*x)*np.cos(np.pi*y)*np.tanh(10000.0*(time-0.01))#np.sin(np.pi * x)* np.cos(np.pi * y) * np.sin(np.pi * 0.5 * time/0.01)
    
    # This function implements the MLS approximation
    def mls_interpolation(self, x, y, values, x_query, k=5, scale_factor=1.5, order=2):
        """Moving Least Squares (MLS) with selectable polynomial order."""

        def weight_function(xi, xj, h):
            """Gaussian weight function for MLS."""
            return np.exp(-np.linalg.norm(xi - xj) ** 2 / (h ** 2))

        def calculate_h(known_points, x_test, k=5, scale_factor=1.5):
            """Estimate local h based on nearest neighbor distances."""
            tree = KDTree(known_points)  
            distances, _ = tree.query(x_test, k+1)  
            return scale_factor * np.mean(distances[1:])  

        def generate_basis(x, y, order):
            """Generate polynomial basis terms dynamically based on order."""
            terms = []
            for i, j in product(range(order + 1), repeat=2):  
                if i + j <= order:  
                    terms.append(x**i * y**j)
            return np.array(terms)
        
        if (len(x) < k):
            k = len(x)

        # Convert data to NumPy arrays
        x, y, values = map(np.array, (x, y, values))

        # Build KDTree for fast nearest neighbor search
        known_points = np.vstack((x, y)).T
        tree = KDTree(known_points)

        # Find k nearest neighbors
        distances, indices = tree.query(x_query, k=k)
        x_k, y_k, values_k = x[indices], y[indices], values[indices]

        # Compute local h dynamically
        h_dynamic = calculate_h(known_points, x_query, k, scale_factor)

        # Generate basis terms for the given order
        basis_size = (order + 1) * (order + 2) // 2  # Compute the number of basis terms
        A = np.zeros((basis_size, basis_size))
        b = np.zeros(basis_size)

        # Compute weighted least squares
        for xi, yi, vi in zip(x_k, y_k, values_k):
            phi = generate_basis(xi, yi, order)  
            w = weight_function(np.array([xi, yi]), np.array(x_query), h_dynamic)
            A += w * np.outer(phi, phi)
            b += w * vi * phi

        # Solve for coefficients
        coeffs = np.linalg.solve(A, b)

        # Compute the interpolated value at x_query
        return coeffs @ generate_basis(x_query[0], x_query[1], order)

    def ComputeErrorAndUpdateOldSolutionField(self, iteration_number):
        # Get the unknown variable for the specific physical problem which is being solved
        unknown_variable_name = self.solver_settings["convection_diffusion_variables"]["unknown_variable"].GetString()
        unknown_variable = KratosMultiphysics.KratosGlobals.GetVariable(unknown_variable_name)

        self.new_solution_field = np.zeros(self._GetSolver().GetComputingModelPart().NumberOfNodes())

        for element in self._GetSolver().GetComputingModelPart().Elements:
            for node in element.GetNodes():
                self.new_solution_field[node.Id - 1] = node.GetSolutionStepValue(unknown_variable)

        error = np.linalg.norm(self.new_solution_field - self.old_solution_field)/np.linalg.norm(self.old_solution_field)

        self.old_solution_field = self.new_solution_field

        return error

    def ApplyDirichletBoundaryConditions(self):
        # Get the model part to apply the dirichlet boundary conditions
        self.model_part = self._GetSolver().GetComputingModelPart()
        self.number_of_nodes_mp = self.model_part.NumberOfNodes()

        # Get the submodel part containing the intersected elements, the active elements and the embedded body
        self.intersected_elements_sub_model_part = self.model_part.GetSubModelPart("intersected_elements")
        self.active_elements_sub_model_part = self.model_part.GetSubModelPart("active_elements")
        self.inactive_elements_sub_model_part = self.model_part.GetSubModelPart("inactive_elements")
        self.embedded_body_boundary_model_part = self.model_part.GetModel().GetModelPart("embedded_body_boundary")

        # Get the unknown variable for the specific physical problem which is being solved
        unknown_variable_name = self.solver_settings["convection_diffusion_variables"]["unknown_variable"].GetString()
        unknown_variable = KratosMultiphysics.KratosGlobals.GetVariable(unknown_variable_name)

        # To find the dirichlet boundary conditions to be applied, we need to solve the problem LHS * phi_dir = RHS

        # Build the LHS and the RHS
        if (self.iteration_number == 0):
            # These contributions can be calculated once per time step as they are not changing during iterations
            self.LHS = self.CalculateLHS()
            self.RHSBoundaryContribution = self.CalculateRHSBoundaryContribution()

        # This contribution needs to be updated in each iteration
        self.RHSTrimmedElementsContribution = self.CalculateRHSTrimmedElementsContribution()

        # Assemble the RHS 
        self.RHS = self.RHSBoundaryContribution + self.RHSTrimmedElementsContribution

        # Solve the linear system
        solution_dir, info = spla.cg(self.LHS, self.RHS)
        #solution_dir = np.linalg.solve(self.LHS, self.RHS)
        
        # Fix the value of the nodes in the "intersected_elements" sub model part
        for element in self.intersected_elements_sub_model_part.Elements:
            for node in element.GetNodes():
                solution_value_to_fix = solution_dir[node.Id - 1]
                node.Fix(unknown_variable)
                node.SetSolutionStepValue(unknown_variable, solution_value_to_fix)

        for element in self.inactive_elements_sub_model_part.Elements:
            for node in element.GetNodes():
                solution_value_to_fix = 0.0
                node.Fix(unknown_variable)
                node.SetSolutionStepValue(unknown_variable, solution_value_to_fix)
    
    def CalculateLHS(self):
        """This method calculates the LHS for the problem that needs to be solved to find the Dirichlet BC's

        LHS = \int_\Gamma N^T * N d\Gamma  + \int_{\Omega_{trim}} (B^T*n)(n^T*B) d{\Omega_{trim}} 
        LHS = LHS_Boundary_Contribution + LHS_Trimmed_Elements_Contribution

        See convection_diffusion_solver.py for more information.
        """

        LHSBoundaryContribution = self.CalculateLHSBoundaryContribution()
        LHSTrimmedElementsContribution = self.CalculateLHSTrimmedElementsContribution()

        LHS = LHSBoundaryContribution + LHSTrimmedElementsContribution

        # Verify complete zero lines and modify the diagonal
        for i in range(LHS.shape[0]):
            if np.all(LHS[i, :] == 0):  # Check if the row is full of ceros
                LHS[i, i] = 1  # Change the value in the diagonal to 1

        return LHS 
    
    
    def CalculateLHSBoundaryContribution(self):
        """ This method calculates the LHS Boundary Contribution:

            LHSBoundaryContribution = \int_\Gamma N^T * N d\Gamma 

        """
        LHSBoundaryContribution = np.zeros((self.number_of_nodes_mp, self.number_of_nodes_mp))

        number_of_nodes_element = 0
        for element in self.intersected_elements_sub_model_part.Elements:
            number_of_nodes_element = element.GetGeometry().PointsNumber()
            break

        for condition in self.embedded_body_boundary_model_part.Conditions:
            condition_geometry = condition.GetGeometry()
            number_of_nodes_condition = condition_geometry.PointsNumber()

            # Get the nodes defining the line segment
            first_node = condition.GetNodes()[0]
            last_node = condition.GetNodes()[1]

            # First and last points of the line 2D condition
            first_point = (first_node.X, first_node.Y)
            last_point = (last_node.X, last_node.Y)

            # Create the line segment from start to end
            line_segment = LineString([first_point, last_point])

            # We need to define which elements intersects this line and we should integrate in each element
            for element in self.intersected_elements_sub_model_part.Elements:
                nodes_x_coordinate = []
                nodes_y_coordinate = []
                for node in element.GetNodes():
                    nodes_x_coordinate.append(node.X)
                    nodes_y_coordinate.append(node.Y)
                    if (len(nodes_x_coordinate) >= 4):
                        break

                # Create a polygon for the element
                points_list = list(zip(nodes_x_coordinate, nodes_y_coordinate))     
                element_polygon = Polygon(points_list)

                # For every intersected element, calculate the local contribution
                if line_segment.intersects(element_polygon):
                    # Find the intersection points
                    intersection_points = line_segment.intersection(element_polygon)

                    if (len(intersection_points.coords) == 1):
                        continue
                    
                    # Initialize the elemental matrix
                    LHSBoundaryContribution_elemental = np.zeros((number_of_nodes_element, number_of_nodes_element)) 
                    
                    # Get the integration method for the line integral
                    weights = [1, 1]
                    integration_points = [-1.0/np.sqrt(3), 1.0/np.sqrt(3.0)]

                    # Line differential
                    ds = np.sqrt((intersection_points.coords[1][0] - intersection_points.coords[0][0])**2 + (intersection_points.coords[1][1] - intersection_points.coords[0][1])**2) / 2.0

                    for gp_index in range(len(integration_points)):

                        # Integration point position in the gauss space, which coincides with the parameter space of the line
                        xi = integration_points[gp_index]

                        # Integration point position in the physical space
                        x_t = (intersection_points.coords[1][0] + intersection_points.coords[0][0]) / 2 + xi * (intersection_points.coords[1][0] - intersection_points.coords[0][0]) / 2
                        y_t = (intersection_points.coords[1][1] + intersection_points.coords[0][1]) / 2 + xi * (intersection_points.coords[1][1]- intersection_points.coords[0][1]) / 2

                        # Integration point position in the parametric space of the rectangular element
                        local_coordinates = element.GetGeometry().PointLocalCoordinates([0,0,0], [x_t, y_t, 0.0])
                        
                        N =  element.GetGeometry().VectorShapeFunctionsValues(np.zeros(2), local_coordinates)
                        
                        LHSBoundaryContribution_elemental += np.outer(N, N) * ds * weights[gp_index]

                    # Assemble the elemental contribution in the global matrix
                    element_connectivities = []
                    for node in element.GetNodes():
                        element_connectivities.append(node.Id)

                    for k in range(len(element_connectivities)):
                        for l in range(len(element_connectivities)):
                            LHSBoundaryContribution[element_connectivities[k] - 1, element_connectivities[l] - 1] += LHSBoundaryContribution_elemental[k, l]

        return LHSBoundaryContribution
    
    def CalculateLHSTrimmedElementsContribution(self):
        """ This method calculates the LHS Trimmed Elements Contribution:

            LHSTrimmedElementsContribution = \int_{\Omega_{trim}} (B^T*n)(n^T*B) d{\Omega_{trim}} 

        """
        LHSTrimmedElementsContribution = np.zeros((self.number_of_nodes_mp, self.number_of_nodes_mp))

        # Iterate over the elements, calculate the local contributions and assemble them to the global matrix
        for element in self.intersected_elements_sub_model_part.Elements:
            element_geometry = element.GetGeometry()
            number_of_nodes_element = element_geometry.PointsNumber()

            # Initialize the elemental matrix
            LHSTrimmedElementsContribution_elemental = np.zeros((number_of_nodes_element, number_of_nodes_element)) 

            # Get the integration method 
            integration_method = element_geometry.GetDefaultIntegrationMethod()
            integration_points, weights = self.DefineIntegrationPoints(integration_method)

            # Vector containing the determinant of the jacobian at each integration point 
            determinant_of_jacobian_vector = element_geometry.DeterminantOfJacobian()

            # Calculate the normal vector of the boundary inside the element
            normal_vector = self.CalculateNormalToBoundary(element)
            normal_vector = -1.0 * np.array(normal_vector[:2])
            
            for gp_index in range(len(integration_points)):
                det_jacobian = determinant_of_jacobian_vector[gp_index]
                
                # Jacobian matrix for each integration point 
                jacobian = element_geometry.Jacobian(gp_index)
                inv_jacobian = np.linalg.inv(jacobian)

                # Matrix containing the gradients in local coordinates. Matrix=[[dN1/dxi dN1/deta], [dN2/dxi dN2/deta]]
                shape_functions_local_gradients = element_geometry.ShapeFunctionsLocalGradients(np.zeros((2, 2)), integration_points[gp_index])
                
                # Shape functions global gradients (with respect to x and y)
                B_matrix = shape_functions_local_gradients @ inv_jacobian
                
                normals_product = np.outer(normal_vector, normal_vector)
                # if (self.model_part.ProcessInfo[KratosMultiphysics.TIME] <= 0.003):
                #     epsilon = 0.001
                # else:
                epsilon = 1.0

                LHSTrimmedElementsContribution_elemental +=  epsilon * (B_matrix @ normals_product @ B_matrix.T) * det_jacobian * weights[gp_index] 
                
            
            # Assemble the elemental contribution in the global matrix
            element_connectivities = []
            for node in element.GetNodes():
                element_connectivities.append(node.Id)

            for k in range(len(element_connectivities)):
                for l in range(len(element_connectivities)):
                    LHSTrimmedElementsContribution[element_connectivities[k] - 1, element_connectivities[l] - 1] += LHSTrimmedElementsContribution_elemental[k, l]
        
        return LHSTrimmedElementsContribution
    
    def CalculateRHSBoundaryContribution(self):
        """ This method calculates the RHS Boundary Contribution:

            RHSBoundaryContribution = \int_\Gamma N * \bar{\phi} d\Gamma 

        """
        RHSBoundaryContribution = np.zeros((self.number_of_nodes_mp, 1))

        number_of_nodes_element = 0
        for element in self.intersected_elements_sub_model_part.Elements:
            number_of_nodes_element = element.GetGeometry().PointsNumber()
            break

        for condition in self.embedded_body_boundary_model_part.Conditions:
            condition_geometry = condition.GetGeometry()
            number_of_nodes_condition = condition_geometry.PointsNumber()

            # Get the nodes defining the line segment
            first_node = condition.GetNodes()[0]
            last_node = condition.GetNodes()[1]

            # First and last points of the line 2D condition
            first_point = (first_node.X, first_node.Y)
            last_point = (last_node.X, last_node.Y)

            # Create the line segment from start to end
            line_segment = LineString([first_point, last_point])

            # We need to define which elements is intersecting this line and we should integrate in each element
            for element in self.intersected_elements_sub_model_part.Elements:
                nodes_x_coordinates = []
                nodes_y_coordinates = []
                for node in element.GetNodes():
                    nodes_x_coordinates.append(node.X)
                    nodes_y_coordinates.append(node.Y)
                    if (len(nodes_x_coordinates) >= 4):
                        break

                # Create a polygon for the element
                points_list = list(zip(nodes_x_coordinates, nodes_y_coordinates))     
                element_polygon = Polygon(points_list)

                # For every intersected element, calculate the local contribution
                if line_segment.intersects(element_polygon):
                    # Find the intersection points
                    intersection_points = line_segment.intersection(element_polygon)

                    if (len(intersection_points.coords) == 1):
                        continue
                    
                    # Initialize the elemental matrix
                    RHSBoundaryContribution_elemental = np.zeros((number_of_nodes_element, 1)) 
                    
                    # Get the integration method for the line integral
                    weights = [1, 1]
                    integration_points = [-1.0/np.sqrt(3), 1.0/np.sqrt(3.0)]

                    # Line differential
                    ds = np.sqrt((intersection_points.coords[1][0] - intersection_points.coords[0][0])**2 + (intersection_points.coords[1][1] - intersection_points.coords[0][1])**2) / 2.0

                    for gp_index in range(len(integration_points)):

                        # Integration point position in the gauss space, which coincides with the parameter space of the line
                        xi = integration_points[gp_index]

                        # Integration point position in the physical space
                        x_t = (intersection_points.coords[1][0] + intersection_points.coords[0][0]) / 2.0 + xi * (intersection_points.coords[1][0] - intersection_points.coords[0][0]) / 2.0
                        y_t = (intersection_points.coords[1][1] + intersection_points.coords[0][1]) / 2.0 + xi * (intersection_points.coords[1][1]- intersection_points.coords[0][1]) / 2.0

                        # Integration point position in the parametric space of the element
                        local_coordinates = element.GetGeometry().PointLocalCoordinates([0,0,0], [x_t, y_t, 0.0])
                        
                        N =  np.array(element.GetGeometry().VectorShapeFunctionsValues(np.zeros(2), local_coordinates))

                        BC_value = self.BC_value(x_t, y_t)
                        
                        RHSBoundaryContribution_elemental += (N * BC_value * ds * weights[gp_index]).reshape(-1, 1)

                    # Assemble the elemental contribution in the global matrix
                    element_connectivities = []
                    for node in element.GetNodes():
                        element_connectivities.append(node.Id)

                    for k in range(len(element_connectivities)):
                        RHSBoundaryContribution[element_connectivities[k] - 1, 0] += RHSBoundaryContribution_elemental[k, 0]
                       
        return RHSBoundaryContribution
    
    def CalculateRHSTrimmedElementsContribution(self):
        """ This method calculates the LHS Trimmed Elements Contribution:

            RHSTrimmedElementsContribution = \int_{\Omega_{trim}} (B^T*n)\big(\nabla \phi|_{extrapolated}^{(k)}*n\big) d{\Omega_{trim}} 

        """
        RHSTrimmedElementsContribution = np.zeros((self.number_of_nodes_mp, 1))

        # Construct the RBF Interpolator
        # Phi values to be used to build the RBF
        solution_derivative_x = []
        solution_derivative_y = []
        x_coordinates = []
        y_coordinates = []

        # Get the unknown variable for the specific physical problem which is being solved
        unknown_variable_name = self.solver_settings["convection_diffusion_variables"]["unknown_variable"].GetString()
        unknown_variable = KratosMultiphysics.KratosGlobals.GetVariable(unknown_variable_name)
        
        for element in self.active_elements_sub_model_part.Elements:
            # Coordinates of the element center in physical space and local space
            element_center = element.GetGeometry().Center()
            element_center_local_space = element.GetGeometry().PointLocalCoordinates([0,0,0], element_center)

            # Get the jacobian and its inverse at the element center
            jacobian = element.GetGeometry().Jacobian(element_center)
            inv_jacobian = np.linalg.inv(jacobian)

            # Get the B matrix in the element center
            shape_functions_local_gradients = element.GetGeometry().ShapeFunctionsLocalGradients(np.zeros((2, 2)), element_center_local_space)
            B_matrix = shape_functions_local_gradients @ inv_jacobian

            solution = []
            is_fixed = False
            for node in element.GetNodes():
                if node.IsFixed(unknown_variable) == True:
                    is_fixed = True
                solution.append(node.GetSolutionStepValue(unknown_variable))

            if is_fixed:
                continue
                
            # Coordinates of the element center
            x_coordinates.append(element_center.X)
            y_coordinates.append(element_center.Y)

            solution_gradients = B_matrix.T @ solution

            solution_derivative_x.append(solution_gradients[0]), solution_derivative_y.append(solution_gradients[1])
            
        # Create RBF interpolator
        #solution_derivative_x_rbf_interpolator = Rbf(x_coordinates, y_coordinates, solution_derivative_x, neigbours=30, function='linear') 
        #solution_derivative_y_rbf_interpolator = Rbf(x_coordinates, y_coordinates, solution_derivative_y, neigbours=30, function='linear')
        # solution_derivative_x_rbf_interpolator = Rbf(x_coordinates, y_coordinates, solution_derivative_x, neigbours=5, function='linear') 
        # solution_derivative_y_rbf_interpolator = Rbf(x_coordinates, y_coordinates, solution_derivative_y, neigbours=5, function='linear')
        # solution_derivative_x_rbf_interpolator = Rbf(x_coordinates, y_coordinates, solution_derivative_x, function='multiquadric') 
        # solution_derivative_y_rbf_interpolator = Rbf(x_coordinates, y_coordinates, solution_derivative_y, function='multiquadric')


        # Iterate over the elements, calculate the local contributions and assemble them to the global matrix
        for element in self.intersected_elements_sub_model_part.Elements:
            element_geometry = element.GetGeometry()
            number_of_nodes_element = element_geometry.PointsNumber()

            # Initialize the elemental matrix
            RHSTrimmedElementsContribution_elemental = np.zeros((number_of_nodes_element, 1)) 

            # Get the integration method 
            integration_method = element_geometry.GetDefaultIntegrationMethod()
            integration_points, weights = self.DefineIntegrationPoints(integration_method)

            # Vector containing the determinant of the jacobian at each integration point 
            determinant_of_jacobian_vector = element_geometry.DeterminantOfJacobian()

            # Calculate the normal vector of the boundary inside the element
            normal_vector = self.CalculateNormalToBoundary(element)
            normal_vector = -1.0 * np.array(normal_vector[:2])

            for gp_index in range(len(integration_points)):
                xi = integration_points[gp_index][0]
                eta = integration_points[gp_index][1]

                # Position of the gauss point in the physical space
                global_coordinates = element_geometry.GlobalCoordinates([xi, eta, 0.0])
                x = global_coordinates[0]
                y = global_coordinates[1]

                # Use the RBF Interpolator to obtain the gradients of the solution
                # derivative_x_interpolated= solution_derivative_x_rbf_interpolator(x, y)
                # derivative_y_interpolated = solution_derivative_y_rbf_interpolator(x, y)
                # Create MLS interpolator
                derivative_x_interpolated = self.mls_interpolation(x_coordinates, y_coordinates, solution_derivative_x, (x, y), 200, 1.5, 2)
                derivative_y_interpolated = self.mls_interpolation(x_coordinates, y_coordinates, solution_derivative_y, (x, y), 200, 1.5, 2)
                # Use the exact gradient
                #derivative_x_interpolated= np.pi*np.cos(np.pi*x)*np.cos(np.pi*y)
                #derivative_y_interpolated = -np.pi*np.sin(np.pi*x)*np.sin(np.pi*y)
                # print(x, y)
                # print(derivative_x_interpolated, derivative_y_interpolated)
                # print(np.pi*np.cos(np.pi*x)*np.cos(np.pi*y), -np.pi*np.sin(np.pi*x)*np.sin(np.pi*y))

                # Obtain the normal gradients in the gauss points
                gauss_point_normal_gradient = derivative_x_interpolated * normal_vector[0] + derivative_y_interpolated * normal_vector[1]
            
                # Get the jacobian and the determinant of the jacobian
                jacobian = element_geometry.Jacobian([xi, eta, 0.0])
                inv_jacobian = np.linalg.inv(jacobian)
                det_jacobian = determinant_of_jacobian_vector[gp_index]

                # Get the B matrix in the gauss point
                shape_functions_local_gradients = element_geometry.ShapeFunctionsLocalGradients(np.zeros((2, 2)), [xi, eta, 0.0])
                B_matrix = shape_functions_local_gradients @ inv_jacobian
                
                epsilon = 1.0
                RHSTrimmedElementsContribution_elemental +=  epsilon * ((B_matrix @ normal_vector) * gauss_point_normal_gradient * det_jacobian * weights[gp_index]).reshape(-1, 1)

            # Assemble the elemental contribution in the global matrix
            element_connectivities = []
            for node in element.GetNodes():
                element_connectivities.append(node.Id)

            for k in range(len(element_connectivities)):
                RHSTrimmedElementsContribution[element_connectivities[k] - 1, 0] += RHSTrimmedElementsContribution_elemental[k, 0]

        return RHSTrimmedElementsContribution


    def DefineIntegrationPoints(self, integration_method):
        if (str(integration_method) == "GeometryData_IntegrationMethod.GI_GAUSS_1"):
            integration_points = [[1.0/3.0, 1.0/3.0, 0.0]]
            weights = [1.0]
            return integration_points, weights 
        elif (str(integration_method) == "GeometryData_IntegrationMethod.GI_GAUSS_2"):
            integration_points = [[-1/np.sqrt(3), -1/np.sqrt(3), 0.0], 
                      [-1/np.sqrt(3), 1/np.sqrt(3), 0.0], 
                      [1/np.sqrt(3), -1/np.sqrt(3), 0.0], 
                      [1/np.sqrt(3), 1/np.sqrt(3), 0.0]]
            weights = [1.0, 1.0, 1.0, 1.0]
            return integration_points, weights 
        elif (str(integration_method) == "GeometryData_IntegrationMethod.GI_GAUSS_3"):
            integration_points = [
                [-np.sqrt(3/5), -np.sqrt(3/5), 0.0], [-np.sqrt(3/5), 0.0, 0.0],
                [-np.sqrt(3/5), np.sqrt(3/5), 0.0], [0.0, -np.sqrt(3/5), 0.0],
                [0.0, 0.0, 0.0], [0.0, np.sqrt(3/5), 0.0],
                [np.sqrt(3/5), -np.sqrt(3/5), 0.0], [np.sqrt(3/5), 0.0, 0.0],
                [np.sqrt(3/5), np.sqrt(3/5), 0.0]
            ]
            weights = [
                (5/9) * (5/9), (5/9) * (8/9), (5/9) * (5/9), (8/9) * (5/9), (8/9) * (8/9), (8/9) * (5/9),
                (5/9) * (5/9), (5/9) * (8/9), (5/9) * (5/9)
            ]
            return integration_points, weights 
        

    def CalculateNormalToBoundary(self, element):
        # Build a polygon using the vertices of the element
        nodes_x_coordinate = []
        nodes_y_coordinate = []

        # Recover the nodes coordinates for an element 
        for node in element.GetNodes():
            nodes_x_coordinate.append(node.X)
            nodes_y_coordinate.append(node.Y)
            if (len(nodes_x_coordinate) >= 4):
                    break

        # Create a polygon for the element
        points_list = list(zip(nodes_x_coordinate, nodes_y_coordinate))     
        element_polygon = Polygon(points_list)

        # Check if a line segment lies inside the element and calculate the normal for this segment 
        line_segments_inside_element = []
        for condition in self.embedded_body_boundary_model_part.Conditions:
            # Start and end points of the line 2D condition
            line_start = (condition.GetNodes()[0].X, condition.GetNodes()[0].Y)
            line_end = (condition.GetNodes()[1].X, condition.GetNodes()[1].Y)

            # Create the line segment from start to end
            line_segment = LineString([line_start, line_end])

            if line_segment.intersects(element_polygon):
                line_segments_inside_element.append(condition)
        
        # For each line segment inside the element, calculate the unit normal
        unit_normals = []
        for condition in line_segments_inside_element:
            unit_normals.append(condition.GetGeometry().UnitNormal())
        
        # Calculate the average normal inside the element
        normals_array = np.array(unit_normals)
        average_normal = np.mean(normals_array, axis=0)
        average_normal_normalized = average_normal / np.linalg.norm(average_normal)
        
        return average_normal_normalized


if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 convection_diffusion_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 convection_diffusion_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = ExtendedGradientMethodConvectionDiffusionAnalysis(model, parameters)
    simulation.Run()
