# Importing the Kratos Library
import KratosMultiphysics

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver
import json
import os

def CreateSolver(model, custom_settings):
    return IgaContactMechanicsSolver(model, custom_settings)

class IgaContactMechanicsSolver(MechanicalSolver):
    """The structural mechanics static solver.

    This class creates the mechanical solvers for static analysis.

    See structural_mechanics_solver.py for more information.
    """
    with open('txt_files/input_data.txt', 'r') as file_in:
        line = file_in.readline().strip()
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
                    
                    
                    
        # REMOVE
        line = file_in.readline().strip()
        line = file_in.readline().strip()
        line = file_in.readline().strip()
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
            file_name = os.path.join(directory, "true_points2.txt")
            if os.path.exists(file_name):
                os.remove(file_name)
            with open(file_name, 'w') as file:
                for condition in skin_model_part2.Conditions :
                    file.write(f"{condition.GetNodes()[0].X} {condition.GetNodes()[0].Y}\n")
                    file.write(f"{condition.GetNodes()[1].X} {condition.GetNodes()[1].Y}\n")
        
        line = file_in.readline().strip()
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
            file_name = os.path.join(directory, "true_points3.txt")
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
    
    
    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "Initializing ...")
        
        # The mechanical solution strategy is created here if it does not already exist.
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        mechanical_solution_strategy = self._GetSolutionStrategy()
        mechanical_solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        mechanical_solution_strategy.Initialize()

        # Printing that inialization is finished
        KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "Finished initialization.")
        
    
    def _CreateSolutionStrategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        
        if analysis_type == "linear":
            mechanical_solution_strategy = self._create_linear_strategy()
        elif analysis_type == "non_linear":
            # Create strategy
            if self.settings["solving_strategy_settings"]["type"].GetString() == "newton_raphson":
                mechanical_solution_strategy = self._create_newton_raphson_strategy()
            elif self.settings["solving_strategy_settings"]["type"].GetString() == "line_search":
                mechanical_solution_strategy = self._create_line_search_strategy()
            elif self.settings["solving_strategy_settings"]["type"].GetString() == "arc_length":
                mechanical_solution_strategy = self._create_arc_length_strategy()

        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\""
            raise Exception(err_msg)
        return mechanical_solution_strategy
    
    def _GetSolutionStrategy(self):
        if not hasattr(self, '_mechanical_solution_strategy'):
            self._mechanical_solution_strategy = self._CreateSolutionStrategy()
        elif not self.settings["builder_and_solver_settings"]["use_block_builder"].GetBool(): # Block builder and solver are unified with MPC and without. In the case of the elimination this could be a problem
            if self.GetComputingModelPart().NumberOfMasterSlaveConstraints() > 0 and not self.mpc_block_builder_initialized:
                self._mechanical_solution_strategy = self._CreateSolutionStrategy()
        return self._mechanical_solution_strategy
        
    
    def InitializeSolutionStep(self):
        
        if self.settings["clear_storage"].GetBool():
            self.Clear()
            self.Initialize() #required after clearing
        self._GetSolutionStrategy().InitializeSolutionStep()
        
    def _GetScheme(self):
        if not hasattr(self, '_solution_scheme'):
            self._solution_scheme = self._CreateScheme()
        return self._solution_scheme

    def _CreateScheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
