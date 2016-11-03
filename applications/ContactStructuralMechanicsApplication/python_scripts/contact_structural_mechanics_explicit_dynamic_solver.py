from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the implicit solver (the explicit one is derived from it)
import structural_mechanics_explicit_dynamic_solver

def CreateSolver(main_model_part, custom_settings):
    return ExplicitMechanicalSolver(main_model_part, custom_settings)

class ExplicitMechanicalSolver(structural_mechanics_explicit_dynamic_solver.ExplicitMechanicalSolver):
    
    ##constructor. the constructor shall only take care of storing the settings 
    ##and the pointer to the main_model part. This is needed since at the point of constructing the 
    ##model part is still not filled and the variables are not yet allocated
    ##
    ##real construction shall be delayed to the function "Initialize" which 
    ##will be called once the model is already filled
    def __init__(self, main_model_part, custom_settings): 
        
        #TODO: shall obtain the computing_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part    
        
        #TODO: remove unnecessary fields for the Explicit solver from the defaults
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "contact_structural_mechanics_explicit_dynamic_solver",
            "echo_level": 0,
            "buffer_size": 2,
            "solution_type": "Dynamic",
            "time_integration_method": "Explicit",
            "scheme_type": "CentralDifferences",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "rotation_dofs": false,
            "pressure_dofs": false,
            "compute_mortar_contact": false,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_step": false,
            "compute_reactions": true,
            "move_mesh_flag": true,
            "clear_storage": false,
            "max_delta_time": 1.0e-5,
            "fraction_delta_time": 0.9,
            "time_step_prediction_level": 0,
            "rayleigh_damping": false,
            "linear_solver_settings":{
                "solver_type": "SuperLU",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
            },
            "problem_domain_sub_model_part_list": ["solid_model_part"],
            "processes_sub_model_part_list": [""]
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        #construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        
        print("Construction of MechanicalSolver finished")

    def AddVariables(self):
        
        structural_mechanics_explicit_dynamic_solver.ExplicitMechanicalSolver.AddVariables(self)
            
        if  self.settings["compute_mortar_contact"].GetBool():
            # Add normal
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
            # Add lagrange multiplier
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER)
            # Add weighted gap
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.WEIGHTED_GAP)
            # Add weighted slip
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.WEIGHTED_SLIP)
            # Add weighted friction
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.WEIGHTED_FRICTION)
            # Add normal augmentation factor
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.NORMAL_AUGMENTATION_FACTOR)
            # Add tangent augmentation factor
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.TANGENT_AUGMENTATION_FACTOR)
            # Auxiliar active
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.AUXILIAR_ACTIVE)
            # Auxiliar slip
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.AUXILIAR_SLIP)
            # Active check factor
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.ACTIVE_CHECK_FACTOR)
        print("::[Mechanical Solver]:: Variables ADDED")
        
    def AddDofs(self):

        structural_mechanics_explicit_dynamic_solver.ExplicitMechanicalSolver.AddDofs(self)
        
        if  self.settings["compute_mortar_contact"].GetBool():
            for node in self.main_model_part.Nodes:
                node.AddDof(KratosMultiphysics.StructuralMechanicsApplication.VECTOR_LAGRANGE_MULTIPLIER_X);
                node.AddDof(KratosMultiphysics.StructuralMechanicsApplication.VECTOR_LAGRANGE_MULTIPLIER_Y);
                node.AddDof(KratosMultiphysics.StructuralMechanicsApplication.VECTOR_LAGRANGE_MULTIPLIER_Z);

        print("::[Mechanical Solver]:: DOF's ADDED")

