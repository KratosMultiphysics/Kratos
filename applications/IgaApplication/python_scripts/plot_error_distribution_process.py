# Kratos imports
import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication 

# External imports 
import numpy as np
import sympy as sp
from colorama import Fore, Style, init
import matplotlib.pyplot as plt

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PlotErrorDistributionProcess(Model, settings["Parameters"])


class PlotErrorDistributionProcess(KM.Process):
    def __init__(self, model, settings):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "background_model_part_name"                 : "please_specify_model_part_name",
                "compute_error_model_part_name"                 : "please_specify_model_part_name",
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
        self.compute_error_model_part_name = self.settings["compute_error_model_part_name"].GetString()
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
        self.compute_error_elements = self.background_model_part.GetSubModelPart(self.compute_error_model_part_name)

        # Initialize lists to store the data
        x_coords = []  # x-coordinates of integration points
        y_coords = []  # y-coordinates of integration points
        error_values = []  # Error values at integration points

        # Iterate over the elements and sum the error 
        for element in self.compute_error_elements.Elements:
            element_geometry = element.GetGeometry()

            nodal_solution = []
            for node in element.GetNodes():
                nodal_solution.append(node.GetSolutionStepValue(self.unknown_variable))

            # Position of the gauss point in the physical space
            global_coordinates = element_geometry.Center()
            x = global_coordinates[0]
            y = global_coordinates[1]

            # Compute the numerical solution at the integration point
            N = np.array(element_geometry.ShapeFunctionsValues())
            numerical_solution = np.dot(nodal_solution, N.T)

            # Compute the local error
            error = (self.analytical_solution_evaluation(x, y) - numerical_solution) 

            # Store the data
            x_coords.append(x)
            y_coords.append(y)
            error_values.append(error)

        # Convert lists to numpy arrays for better handling
        x_coords = np.array(x_coords)
        y_coords = np.array(y_coords)
        error_values = np.array(error_values)

        # Create the scatter plot
        plt.figure(figsize=(8, 6))
        scatter = plt.scatter(x_coords, y_coords, c=error_values, cmap="viridis", marker='o')

        # Add color bar to interpret the error values
        plt.colorbar(scatter, label="Error Magnitude")

        # Add labels and title
        plt.xlabel("x-coordinate")
        plt.ylabel("y-coordinate")
        plt.title("Error Distribution Plot")

        # Show the plot
        plt.grid(True, linestyle="--", alpha=0.5)

        # Adding the legend
        plt.legend()
        # Save and show the plot
        plt.savefig('error_distribution_plot.png')
        plt.show()

    def analytical_solution_evaluation(self, x, y):
        return eval(self.analytical_solution)

        

        
        
        
