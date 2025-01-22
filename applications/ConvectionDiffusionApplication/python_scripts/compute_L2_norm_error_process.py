# Kratos imports
import KratosMultiphysics as KM
import KratosMultiphysics.ConvectionDiffusionApplication as CDA

# External imports 
import numpy as np
import sympy as sp
from colorama import Fore, Style, init
import matplotlib.pyplot as plt

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeL2NormErrorProcess(Model, settings["Parameters"])


class ComputeL2NormErrorProcess(KM.Process):
    def __init__(self, model, settings):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "background_model_part_name"                 : "please_specify_model_part_name",
                "analytical_solution"                       : "",
                "unknown_variable_name"                       : ""                           
            }
            """
            )
        
        self.settings = settings

        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model = model


        # Retrieve from the input file the model part names 
        self.background_model_part_name = self.settings["background_model_part_name"].GetString()
        self.analytical_solution = self.settings["analytical_solution"].GetString()
        self.unknown_variable_name = self.settings["unknown_variable_name"].GetString()
        self.unknown_variable = KM.KratosGlobals.GetVariable(self.unknown_variable_name)
        
        # Get the model parts 
        self.background_model_part = self.model.GetModelPart(self.background_model_part_name)


    def Check(self):
        pass

    def ExecuteFinalize(self):
        L2_norm_error = 0.0

        # Get the active elements sub model part
        self.active_elements = self.background_model_part.GetSubModelPart("active_elements")
        # self.active_elements = self.background_model_part


        # Integration points and weights for the numerical integration - Quadrilateral element
        # integration_points = [[-1/np.sqrt(3), -1/np.sqrt(3), 0.0], 
        #               [-1/np.sqrt(3), 1/np.sqrt(3), 0.0], 
        #               [1/np.sqrt(3), -1/np.sqrt(3), 0.0], 
        #               [1/np.sqrt(3), 1/np.sqrt(3), 0.0]]
        # weights = [1.0, 1.0, 1.0, 1.0]

        # Integration points and weights for the numerical integration - Linear triangle element
        integration_points = [[1/2, 1/2, 0.0], [1/2, 0.0, 0.0], [0.0, 1/2, 0.0]]
        weights = [1/3, 1/3, 1/3]

        circumradius = 0.0

        # Iterate over the elements and sum the error 
        for element in self.active_elements.Elements:
            element_geometry = element.GetGeometry()

            nodal_solution = []
            for node in element.GetNodes():
                nodal_solution.append(node.GetSolutionStepValue(self.unknown_variable))

            for gp_index in range(len(integration_points)):
                xi = integration_points[gp_index][0]
                eta = integration_points[gp_index][1]

                # Position of the gauss point in the physical space
                global_coordinates = element_geometry.GlobalCoordinates([xi, eta, 0.0])
                x = global_coordinates[0]
                y = global_coordinates[1]
                t = self.background_model_part.ProcessInfo[KM.TIME]

                # Calculate the jacobian and its determinant at the gauss point
                jacobian = element_geometry.Jacobian([xi, eta, 0.0])
                det_jacobian = np.linalg.det(jacobian)

                # Compute the numerical solution at the integration point
                N = element_geometry.VectorShapeFunctionsValues(np.zeros(2), [xi, eta, 0.0])
                numerical_solution = np.dot(nodal_solution, N)

                # Compute the local error
                error = (self.analytical_solution_evaluation(x, y, t) - numerical_solution) ** 2

                # Contribution to the L2 error 
                L2_norm_error += error * det_jacobian * weights[gp_index]

                circumradius += element_geometry.Circumradius()
        
        print(circumradius/self.active_elements.NumberOfElements())
        print(f"{Fore.GREEN}L2-Norm Error: {np.sqrt(L2_norm_error)}{Style.RESET_ALL}")
        # Data (LINEAR TRIANGLE)
        sizes_emb = [0.0964, 0.0707, 0.0424, 0.0316, 0.028, 0.023, 0.022, 0.0202]
        errors_emb = [0.00614, 0.0029, 0.00094, 0.00050, 0.000425, 0.00027, 0.000219, 0.000189]

        # Data(QUADRATIC QUADRILATERAL)
        # sizes_emb = [0.10, 0.067, 0.05, 0.033, 0.02, 0.014, 0.011, 0.0091]
        # errors_emb = [0.025, 0.00765, 0.0011, 0.000805, 0.00020, 0.0001, 6.01e-5, 1.3e-5]

        # sizes_bf = [0.07543, 0.0628, 0.0535, 0.0332, 0.02915]
        # errors_bf = [8.81e-5, 8.435e-5, 8.231e-5, 3.718e-6, 2.45e-6]

        # Create the plot
        plt.figure(figsize=(6, 4))

        # Set log-log scale
        plt.loglog(sizes_emb, errors_emb, marker='o', linestyle='-', color='b',label="Embedded error")
        # plt.loglog(sizes_bf, errors_bf, marker='x', linestyle='-', color='r', label="body-fitted error")

        # Define the C constant for the slope 3 line
        C = errors_emb[0] / sizes_emb[0]**2 + 1.0  # Usamos el primer punto para estimar C

        # Generate x values for the line
        x_values = np.logspace(np.log10(min(sizes_emb)), np.log10(max(sizes_emb)), 100) # Valores de x en escala logarítmica

        # Calculate the y values
        y_values = C * x_values**2

        # Plot the line
        plt.loglog(x_values, y_values, linestyle='--', color='r', label="Slope 2")

        # Do a lineal fitting through the points
        log_h = np.log(sizes_emb)
        log_error = np.log(errors_emb)
        slope, intercept = np.polyfit(log_h, log_error, 1)  # Retorna la pendiente y el intercepto
        plt.loglog(sizes_emb, np.exp(intercept) * sizes_emb**slope, '--', label=f'Lineal Fitting (slope = {slope:.2f})')

        # Labels and title
        plt.xlabel('Size')
        plt.ylabel('Error')
        plt.title('Log-Log Plot of Size vs Error')

        # Show the grid
        plt.grid(True, which="both", ls="--")
        # Add the legend
        plt.legend()
        # Save and show the plot
        plt.savefig('L2_Norm_error_plot.png')
        plt.show()

    def analytical_solution_evaluation(self, x, y, t):
        return eval(self.analytical_solution)

        

        
        
        
