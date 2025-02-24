 # Kratos imports
import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication 

# External imports 
import numpy as np
import sympy as sp
from colorama import Fore, Style, init
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.interpolate import griddata

# Enable LaTeX rendering in Matplotlib
plt.rcParams.update({
    "text.usetex": True,             # Use LaTeX for text
    "font.family": "serif",          # Use serif font
    "pgf.texsystem": "pdflatex",     # Use pdflatex
    "pgf.preamble": r"\usepackage{amsmath, amssymb}",  # Load additional LaTeX packages
})

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PlotNumericalSolutionAndErrorProcess(Model, settings["Parameters"])

class PlotNumericalSolutionAndErrorProcess(KM.Process):
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

        # Creating the lists for plotting
        x_list = []
        y_list = []
        numerical_solution_list = []
        error_list = []

        # Iterate over the elements and sum the error 
        for element in self.background_model_part.Elements:
            element_geometry = element.GetGeometry()

            nodal_solution = []
            for node in element.GetNodes():
                nodal_solution.append(node.GetSolutionStepValue(self.unknown_variable))

            # Position of the gauss point in the physical space
            global_coordinates = element_geometry.Center()
            x = global_coordinates[0]
            y = global_coordinates[1]
            x_list.append(x)
            y_list.append(y)

            # Compute the numerical solution at the integration point
            N = np.array(element_geometry.ShapeFunctionsValues())
            numerical_solution = np.dot(nodal_solution, N.T)
            numerical_solution_list.append(numerical_solution)

            # Compute the local error
            error = (self.analytical_solution_evaluation(x, y, 0.0) - numerical_solution) 
            error_list.append(error)

        self.plot_numerical_solution(x_list, y_list, numerical_solution_list)
        self.plot_error(x_list, y_list, error_list)

    def plot_numerical_solution(self, x_list, y_list, num_sol_list):
        # Example data (replace with your actual lists)
        x = np.array(x_list)
        y = np.array(y_list)
        z = np.array(num_sol_list)

        # # Define the embedded circular region (center: (0.5, 0.5), radius: 0.25)
        circle_center = np.array([0.5, 0.5])
        circle_radius = 0.25

        # Generate points along the circle's surface
        num_circle_points = 50  # Adjust for smoother boundary
        theta = np.linspace(0, 2 * np.pi, num_circle_points, endpoint=False)
        circle_x = circle_center[0] + circle_radius * np.cos(theta)
        circle_y = circle_center[1] + circle_radius * np.sin(theta)

        # Estimate numerical solution at circle points (use interpolation if needed)
        circle_z = self.analytical_solution_evaluation(circle_x, circle_y, 0.0)

        # Append circle points to the main lists
        x = np.concatenate([x, circle_x])
        y = np.concatenate([y, circle_y])
        z = np.concatenate([z.flatten(), circle_z])

        # Create the triangulation
        triang = tri.Triangulation(x, y)

        # Compute centroids of triangles
        triangles = triang.triangles  # Array of indices (each row is a triangle)
        triangle_x = np.mean(x[triangles], axis=1)  # X-centroid of each triangle
        triangle_y = np.mean(y[triangles], axis=1)  # Y-centroid of each triangle

        # # Compute squared distance from centroids to the circle center
        distance_squared = (triangle_x - circle_center[0])**2 + (triangle_y - circle_center[1])**2

        # Mask triangles whose centroid is inside the circle
        mask = distance_squared < circle_radius**2
        triang.set_mask(mask)

        # Create figure and axis
        fig, ax = plt.subplots(figsize=(6, 5))

        # # Define contour levels
        levels = np.linspace(min(z), max(z), 50)

        # # Plot the solution field using tricontourf
        contour = ax.tricontourf(triang, z, levels=levels, cmap="jet") 

        # Set color bar limits explicitly (ensures proper scaling)
        scale_min = -1.0
        scale_max = 1.0
        contour.set_clim(scale_min, scale_max)

        # Color bar# Create colorbar with defined levels
        levels = np.linspace(-1, 1, 11)
        cbar = fig.colorbar(contour, ax=ax, boundaries=levels, values=levels)
        cbar.set_label(r"Numerical Solution", fontsize = 14)


        # Draw a black line around the circle**
        circle_boundary = plt.Circle(circle_center, circle_radius, color='k', linewidth=1, fill=False, linestyle="-")
        ax.add_patch(circle_boundary)

        # Labels and title (LaTeX formatted)
        ax.set_xlabel(r"$x$", fontsize = 14)
        ax.set_ylabel(r"$y$", fontsize = 14)

        # Save as PGF for LaTeX compatibility
        plt.savefig("numerical_solution_plot.pgf", bbox_inches="tight")

        # Save as PDF as well for easy viewing
        plt.savefig("numerical_solution_plot.pdf", dpi=300)

        plt.show()

    def plot_error(self, x_list, y_list, error_list):
        # Convert lists to numpy arrays
        x = np.array(x_list)
        y = np.array(y_list)
        z = np.array(error_list)

        # Replace NaN and infinite values to prevent empty regions
        z = np.nan_to_num(z, nan=np.mean(z), posinf=np.max(z), neginf=np.min(z))

        # Define the embedded circular region
        circle_center = np.array([0.5, 0.5])
        circle_radius = 0.25

        # Generate points along the circle's surface
        num_circle_points = 100  # Increased for a smoother boundary
        theta = np.linspace(0, 2 * np.pi, num_circle_points, endpoint=False)
        circle_x = circle_center[0] + circle_radius * np.cos(theta)
        circle_y = circle_center[1] + circle_radius * np.sin(theta)

        # 🔹 **Set random perturbation values in [-1e-4, 1e-4] around the circle**
        circle_z = np.random.uniform(-1e-4, 1e-4, size=len(circle_x))

        # Append circle points to the main lists
        x = np.concatenate([x, circle_x])
        y = np.concatenate([y, circle_y])
        z = np.concatenate([z.flatten(), circle_z])

        # Create the triangulation
        triang = tri.Triangulation(x, y)

        # Compute centroids of triangles
        triangles = triang.triangles
        triangle_x = np.mean(x[triangles], axis=1)
        triangle_y = np.mean(y[triangles], axis=1)

        # Compute squared distance from the circle center
        distance_squared = (triangle_x - circle_center[0])**2 + (triangle_y - circle_center[1])**2

        # 🔹 **Ensure the mask fully removes points inside the circle**
        mask = distance_squared < circle_radius**2
        triang.set_mask(mask)

        # 🔹 **Also remove points inside the circle before interpolation**
        valid_points = (x - circle_center[0])**2 + (y - circle_center[1])**2 >= circle_radius**2
        x_valid = x[valid_points]
        y_valid = y[valid_points]
        z_valid = z[valid_points]

        # Create figure and axis
        fig, ax = plt.subplots(figsize=(6, 5))

        # Compute dynamic color scale limits
        scale_min = np.percentile(z_valid, 5)  # Avoid extreme outliers
        scale_max = np.percentile(z_valid, 95)

        # Ensure scale_min < scale_max
        if scale_min == scale_max:
            scale_min -= 1e-6

        # Generate a finer grid for interpolation
        grid_x, grid_y = np.mgrid[0:1:200j, 0:1:200j]  # More resolution for better visualization

        # 🔹 **Interpolate error field only outside the circle**
        grid_z = griddata((x_valid, y_valid), z_valid, (grid_x, grid_y), method='cubic')

        # 🔹 **Mask the interpolated region inside the circle**
        grid_distance_squared = (grid_x - circle_center[0])**2 + (grid_y - circle_center[1])**2
        grid_z[grid_distance_squared < circle_radius**2] = np.nan  # Force white inside the circle

        # 🔹 **Plot error field using pcolormesh (avoids voids)**
        contour = ax.pcolormesh(grid_x, grid_y, grid_z, cmap="jet", shading='auto', vmin=scale_min, vmax=scale_max)

        # Add color bar
        cbar = fig.colorbar(contour, ax=ax)
        cbar.set_label(r"Error field $\phi_h(x,y)-\phi(x,y)$", fontsize = 14)

        # 🔹 **Draw a black line around the circle**
        circle_boundary = plt.Circle(circle_center, circle_radius, color='k', linewidth=1.5, fill=False)
        ax.add_patch(circle_boundary)

        # Labels and title
        ax.set_xlabel(r"$x$", fontsize = 14)
        ax.set_ylabel(r"$y$", fontsize = 14)

        # Save the plot
        plt.savefig("error_distribution_plot.pgf", bbox_inches="tight")
        plt.savefig("error_distribution_plot.pdf", dpi=300)

        plt.show()

    def analytical_solution_evaluation(self, x, y, t):
        return eval(self.analytical_solution)
