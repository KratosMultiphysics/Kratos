# Importing the Kratos Library
import KratosMultiphysics

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

import KratosMultiphysics.IgaApplication as IgaApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import json
import os

# Other imports
from KratosMultiphysics.StructuralMechanicsApplication import IGA_convergence_criteria_factory 

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
                    
    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "solver_type" : "mechanical_solver",
            "model_part_name" : "",
            "computing_sub_model_part_name" : "",
            "domain_size" : -1,
            "echo_level": 0,
            "buffer_size": 2,
            "analysis_type": "non_linear",
            "model_import_settings": {
                "input_type": "mdpa"
            },
            "material_import_settings" :{
                "materials_filename": ""
            },
            "time_stepping" : { },
            "volumetric_strain_dofs": false,
            "rotation_dofs": false,
            "pressure_dofs": false,
            "displacement_control": false,
            "reform_dofs_at_each_step": false,
            "use_old_stiffness_in_first_iteration": false,
            "compute_reactions": true,
            "solving_strategy_settings": {
                "type" : "newton_raphson",
                "advanced_settings" : { }
            },
            "builder_and_solver_settings" : {
                "use_block_builder" : true,
                "use_lagrange_BS"   : false,
                "advanced_settings" : { }
            },
            "clear_storage": false,
            "move_mesh_flag": true,
            "multi_point_constraints_used": true,
            "convergence_criterion": "residual_criterion",
            "scheme_type": "contact",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 10,
            "linear_solver_settings": { },
            "auxiliary_variables_list" : [],
            "auxiliary_dofs_list" : [],
            "auxiliary_reaction_list" : [],
            
            "iga_contact_parameters" : []
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults
                
    
    
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
        """Create the scheme for the scipy solver.

        The scheme determines the mass and stiffness matrices
        """
        scheme_type = self.settings["scheme_type"].GetString()
        
        # solution_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        if scheme_type == "iga_contact":
            solution_scheme = IgaApplication.IgaContactScheme()
            # solution_scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
        else: # here e.g. a stability scheme could be added
            err_msg =  "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"dynamic\""
            # raise Exception(err_msg)

        return solution_scheme
        # return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

    def _create_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._GetScheme()
        mechanical_convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        
        # strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(computing_model_part,
        #                                                              mechanical_scheme,
        #                                                              mechanical_convergence_criterion,
        #                                                              builder_and_solver,
        #                                                              self.settings["max_iteration"].GetInt(),
        #                                                              self.settings["compute_reactions"].GetBool(),
        #                                                              self.settings["reform_dofs_at_each_step"].GetBool(),
        #                                                              self.settings["move_mesh_flag"].GetBool())
        
        strategy = IgaApplication.IgaResidualBasedNewtonRaphsonContactStrategy(computing_model_part,
                                                                            mechanical_scheme,
                                                                            mechanical_convergence_criterion,
                                                                            builder_and_solver,
                                                                            self.settings["iga_contact_parameters"],
                                                                            self.settings["max_iteration"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool())
        strategy.SetUseOldStiffnessInFirstIterationFlag(self.settings["use_old_stiffness_in_first_iteration"].GetBool())
        return strategy
    
    
    def _CreateConvergenceCriterion(self):
        convergence_criterion = IGA_convergence_criteria_factory.convergence_criterion(self._get_convergence_criterion_settings())
        return convergence_criterion.mechanical_convergence_criterion