# Importing the Kratos Library
import KratosMultiphysics

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver
import json
import os

def CreateSolver(model, custom_settings):
    return StaticMechanicalSolver(model, custom_settings)

class StaticMechanicalSolver(MechanicalSolver):
    """The structural mechanics static solver.

    This class creates the mechanical solvers for static analysis.

    See structural_mechanics_solver.py for more information.
    """
    with open('txt_files/input_data.txt', 'r') as file:
        line = file.readline().strip()
    name_mdpa_true_boundary = line
    # name_mdpa_true_boundary = "mdpa_files/Weird_shape3" 
    file_mdpa_exists = os.path.isfile(os.path.join(name_mdpa_true_boundary + ".mdpa"))
    
    if (file_mdpa_exists) :
        print('mdpa file of embedded boundary exists!')
        current_model = KratosMultiphysics.Model()
        skin_model_part2 = current_model.CreateModelPart("skin_model_part2")
        skin_model_part2.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        KratosMultiphysics.ModelPartIO(name_mdpa_true_boundary).ReadModelPart(skin_model_part2)
        
        # Write all the points of the skin boundary in an external file
        directory = "txt_files"
        file_name = os.path.join(directory, "true_points.txt")
        if os.path.exists(file_name):
            os.remove(file_name)
        with open(file_name, 'w') as file:
            for condition in skin_model_part2.Conditions :
                file.write(f"{condition.GetNodes()[0].X} {condition.GetNodes()[0].Y}\n")
                file.write(f"{condition.GetNodes()[1].X} {condition.GetNodes()[1].Y}\n")
                
    
    
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[StaticMechanicalSolver]:: ", "Construction finished")

    def _CreateScheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
