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
        self.compute_error_elements = self.background_model_part.GetSubModelPart("compute_error")
        # self.compute_error_elements = self.model.GetModelPart("ThermalModelPart")

        # Iterate over the elements and sum the error 
        for element in self.compute_error_elements.Elements:
            element_geometry = element.GetGeometry()

            # Get the integration weight
            weight = element.GetValue(KratosMultiphysics.IgaApplication.INTEGRATION_WEIGHTS)

            nodal_solution = []
            skip_element = False
            for node in element.GetNodes():
                # Check if the node is fixed
                # if node.IsFixed(self.unknown_variable):
                #     skip_element = True
                #     break 
                nodal_solution.append(node.GetSolutionStepValue(self.unknown_variable))

            if skip_element:
                continue  # Skip to the next element in the loop


            # Position of the gauss point in the physical space
            global_coordinates = element_geometry.Center()
            x = global_coordinates[0]
            y = global_coordinates[1]

            # Calculate the determinant of the jacobian 
            determinant_of_jacobian = element_geometry.DeterminantOfJacobian()[0]

            # Compute the numerical solution at the integration point
            N = np.array(element_geometry.ShapeFunctionsValues())
            numerical_solution = np.dot(nodal_solution, N.T)

            # Compute the local error
            error = (self.analytical_solution_evaluation(x, y) - numerical_solution) ** 2

            # Contribution to the L2 error 
            L2_norm_error += error * determinant_of_jacobian * weight[0]

                # circumradius += element_geometry.Circumradius()
        
        # print(circumradius/self.compute_error_elements.NumberOfElements())
        print(f"{Fore.GREEN}L2-Norm Error: {np.sqrt(L2_norm_error)}{Style.RESET_ALL}")
        # # Linear errors - MLS
        # sizes_emb_EGM = [ 0.0909, 0.0476, 0.0322, 0.0244, 0.0196, 0.0164, 0.0141, 0.0123, 0.0099]
        # errors_emb_EGM = [0.0244, 0.00157, 0.00034, 0.000103, 4.037e-5, 2.467e-5, 1.729e-5, 1.299e-5, 8.165e-6]
        # # IBRA Errors - Linear
        # sizes_emb_IBRA = [0.1, 0.05, 0.033, 0.025, 0.0227, 0.0163]
        # errors_emb_IBRA = [0.00048, 3.931e-5, 7.533e-6, 2.773e-6, 1.856e-6, 4.88e-7]
        # # Exact Gradient errors - EGM
        # sizes_emb_EGM_exact_grad = [ 0.0909, 0.0476, 0.0322, 0.0244, 0.0178, 0.0163, 0.0131, 0.0109, 0.0099]
        # errors_emb_EGM_exact_grad = [0.000165, 2.284e-5, 5.435e-6, 1.953e-6, 5.911e-7, 4.256e-7, 1.811e-7, 8.853e-8, 5.941e-8]
       
        # # Quadratic errors (50 points) - Linear MLS
        # sizes_emb_EGM = [0.0416, 0.03125, 0.025, 0.0208, 0.0156, 0.0125]
        # errors_emb_EGM = [0.0023, 0.00059, 0.00034, 0.00018, 9.81e-5, 6.198e-5]
        # # Quadratic errors (50 points) - Quadratic MLS
        # sizes_emb_EGM = [0.0625, 0.0416, 0.03125, 0.025, 0.0208, 0.0139, 0.0096]
        # errors_emb_EGM = [0.0055, 0.00043, 0.000138, 8.58e-5, 5.106e-5, 2.18e-5, 6.244e-6]
        # # # IBRA Errors
        # sizes_emb_IBRA = [0.166, 0.0909, 0.0476, 0.0322, 0.0243, 0.0217, 0.0163]
        # errors_emb_IBRA = [2.861e-5, 1.561e-6, 1.88e-7, 3.696e-8, 1.376e-8, 1.216e-8, 4.593e-9]
        # # Quadratic errors (50 points) -- exact gradient
        # sizes_emb_EGM_exact_grad = [0.125, 0.0625, 0.0417, 0.03125, 0.025, 0.0208, 0.0156, 0.0125]
        # errors_emb_EGM_exact_grad = [3.913e-7, 7.667e-9, 7.151e-10, 1.287e-10, 3.441e-11, 1.19e-11, 2.117e-12, 5.551e-13]


        #cubic errors - MLS (linear MLS)
        sizes_emb_EGM = [0.0625, 0.03125, 0.0238, 0.0192, 0.0161, 0.0139]
        errors_emb_EGM_old = [0.005, 0.00028, 0.00012, 5.188e-5, 2.037e-5, 1.366e-5] 
        errors_emb_EGM = [x - 1e-5 for x in errors_emb_EGM_old]
        #cubic errors - MLS (quadratic MLS)
        # sizes_emb_EGM = [0.0625, 0.0454, 0.03125, 0.0238, 0.0192, 0.0139, 0.0122, 0.0108]
        # errors_emb_EGM = [0.031, 0.0045, 0.00056, 0.00032, 0.00013, 4.517e-5, 2.8e-5, 1.90e-5]
        # # IBRA Errors
        sizes_emb_IBRA = [0.0909, 0.0476, 0.0285, 0.0217, 0.0166]
        errors_emb_IBRA = [2.204e-6, 2.057e-7, 2.409e-8, 7.06e-9, 4.491e-9]
        # Quadratic errors (50 points) -- exact gradient
        sizes_emb_EGM_exact_grad = [0.1, 0.0625, 0.05, 0.0384, 0.03125, 0.0238, 0.0192]
        errors_emb_EGM_exact_grad = [3.21e-9, 7.17e-11, 1.617e-11, 7e-13, 9.49e-14, 3.3026e-14, 1.678e-15]



        # Create the figure
        plt.figure(figsize=(8, 11))

        # Set log-log scale
        plt.loglog(sizes_emb_EGM, errors_emb_EGM, marker='o', linestyle='-', color='b',label="New Method Error - MLS Gradient Approximation")
        plt.loglog(sizes_emb_IBRA, errors_emb_IBRA, marker='x', linestyle='-', color='r', label="IBRA Error")
        plt.loglog(sizes_emb_EGM_exact_grad, errors_emb_EGM_exact_grad, marker='o', linestyle='-', color='g',label="New Method Error - Exact Gradient")


        # Define the c constant for a slope of 2
        C = errors_emb_EGM[0] / sizes_emb_EGM[0]**4 - 170.0
        # Generate x values
        x_values = np.logspace(np.log10(min(sizes_emb_EGM)), np.log10(max(sizes_emb_EGM)), 100)
        # Calculate the y values
        y_values = C * x_values**4
        # Plot the slope 3 line
        plt.loglog(x_values, y_values, linestyle='--', color='m', label="Slope 4")

        # Do a linear fitting through the error points
        log_h = np.log(sizes_emb_EGM)
        log_error = np.log(errors_emb_EGM)
        slope, intercept = np.polyfit(log_h, log_error, 1)  
        plt.loglog(sizes_emb_EGM, np.exp(intercept) * sizes_emb_EGM**slope, '--', label=f'Lineal Fitting (slope = {slope:.2f})')

        # Add more numbers on the x-axis
        ax = plt.gca()
        x_ticks = np.logspace(np.log10(0.01), np.log10(0.10), 7)  # Add 15 evenly spaced ticks
        ax.set_xticks(x_ticks)
        ax.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.4f}'))  # Format to 4 decimal places

        # Labels and title
        plt.xlabel('Size')
        plt.ylabel('Error')
        plt.title('Comparative Study: Log-Log Plot Size vs Error')

        # Show the grid
        plt.grid(True, which="both", ls="--")
        # Adding the legend
        plt.legend()
        # Save and show the plot
        plt.savefig('L2_Norm_error_plot.png')
        plt.show()

    def analytical_solution_evaluation(self, x, y):
        return eval(self.analytical_solution)

        

        
        
        
