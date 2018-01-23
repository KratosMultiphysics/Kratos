from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the implicit solver (the explicit one is derived from it)
import structural_mechanics_implicit_dynamic_solver

def CreateSolver(main_model_part, custom_settings):
    return ImplicitMechanicalSolver(main_model_part, custom_settings)

class ImplicitMechanicalSolver(structural_mechanics_implicit_dynamic_solver.ImplicitMechanicalSolver):
    """The structural mechanics contact implicit dynamic solver.

    This class creates the mechanical solvers for contact implicit dynamic analysis.
    It currently supports Newmark, Bossak and dynamic relaxation schemes.

    Public member variables:
    dynamic_settings -- settings for the implicit dynamic solvers.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings): 
        
        self.main_model_part = main_model_part    
        
        ##settings string in json format
        contact_settings = KratosMultiphysics.Parameters("""
        {
            "contact_settings" :
            {
                "mortar_type"                            : "",
                "condn_convergence_criterion"            : false,
                "fancy_convergence_criterion"            : true,
                "print_convergence_criterion"            : false,
                "ensure_contact"                         : false,
                "gidio_debug"                            : false,
                "adaptative_strategy"                    : false,
                "split_factor"                           : 10.0,
                "max_number_splits"                      : 3,
                "contact_displacement_relative_tolerance": 1.0e-4,
                "contact_displacement_absolute_tolerance": 1.0e-9,
                "contact_residual_relative_tolerance"    : 1.0e-4,
                "contact_residual_absolute_tolerance"    : 1.0e-9
            }
        }
        """)
        
        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.validate_and_transfer_matching_settings(self.settings, contact_settings)
        self.contact_settings = contact_settings["contact_settings"]

        # Construct the base solver.
        super().__init__(self.main_model_part, self.settings)
        
        # Setting default configurations true by default
        if (self.settings["clear_storage"].GetBool() == False):
            print("WARNING:: Storage must be cleared each step. Switching to True")
            self.settings["clear_storage"].SetBool(True)
        if (self.settings["reform_dofs_at_each_step"].GetBool() == False):
            print("WARNING:: DoF must be reformed each time step. Switching to True")
            self.settings["reform_dofs_at_each_step"].SetBool(True)
        if (self.settings["block_builder"].GetBool() == False):
            print("WARNING:: EliminationBuilderAndSolver can not used with the current implementation. Switching to BlockBuilderAndSolver")
            self.settings["block_builder"].SetBool(True)
        
        # Setting echo level
        self.echo_level =  self.settings["echo_level"].GetInt()
    
        # Initialize the processes list
        self.processes_list = None
        
        print("Construction of ContactMechanicalSolver finished")

    def AddVariables(self):
        
        super().AddVariables()
            
        if  self.contact_settings["mortar_type"].GetString() != "":
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)  # Add normal
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H) # Add nodal size variable
            if  self.contact_settings["mortar_type"].GetString() == "ALMContactFrictionless":
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_CONTACT_STRESS)                        # Add normal contact stress
                self.main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.WEIGHTED_GAP)              # Add normal contact gap
            elif self.contact_settings["mortar_type"].GetString() == "ALMContactFrictional": 
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER)                   # Add normal contact stress 
                self.main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.WEIGHTED_GAP)              # Add normal contact gap 
                self.main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.WEIGHTED_SLIP)             # Add normal contact gap 
            elif  self.contact_settings["mortar_type"].GetString() == "ScalarMeshTying":
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SCALAR_LAGRANGE_MULTIPLIER)                   # Add scalar LM
                self.main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.WEIGHTED_SCALAR_RESIDUAL)  # Add scalar LM residual
            elif  self.contact_settings["mortar_type"].GetString() == "ComponentsMeshTying":
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER)                   # Add vector LM
                self.main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.WEIGHTED_VECTOR_RESIDUAL)  # Add vector LM residual
   
        print("::[Contact Mechanical Solver]:: Variables ADDED")
        
    def AddDofs(self):

        super().AddDofs()
        
        if (self.contact_settings["mortar_type"].GetString() == "ALMContactFrictionless"):
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.NORMAL_CONTACT_STRESS, ContactStructuralMechanicsApplication.WEIGHTED_GAP, self.main_model_part)
        elif (self.contact_settings["mortar_type"].GetString() == "ALMContactFrictional"): 
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_X, self.main_model_part) 
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_Y, self.main_model_part) 
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_Z, self.main_model_part) 
        elif (self.contact_settings["mortar_type"].GetString() == "ScalarMeshTying"):
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.SCALAR_LAGRANGE_MULTIPLIER,ContactStructuralMechanicsApplication.WEIGHTED_SCALAR_RESIDUAL, self.main_model_part)
        elif (self.contact_settings["mortar_type"].GetString() == "ComponentsMeshTying"):
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_X, ContactStructuralMechanicsApplication.WEIGHTED_VECTOR_RESIDUAL_X, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_Y, ContactStructuralMechanicsApplication.WEIGHTED_VECTOR_RESIDUAL_Y, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_Z, ContactStructuralMechanicsApplication.WEIGHTED_VECTOR_RESIDUAL_Z, self.main_model_part)

        print("::[Contact Mechanical Solver]:: DOF's ADDED")
    
    def Initialize(self):
        super().Initialize() # The mechanical solver is created here.
    
    def AddProcessesList(self, processes_list):
        self.processes_list = ContactStructuralMechanicsApplication.ProcessFactoryUtility(processes_list)
    
    def _create_convergence_criterion(self):
        # Create an auxiliary Kratos parameters object to store the convergence settings.
        conv_params = KratosMultiphysics.Parameters("{}")
        conv_params.AddValue("convergence_criterion",self.settings["convergence_criterion"])
        conv_params.AddValue("rotation_dofs",self.settings["rotation_dofs"])
        conv_params.AddValue("echo_level",self.settings["echo_level"])
        conv_params.AddValue("displacement_relative_tolerance",self.settings["displacement_relative_tolerance"])
        conv_params.AddValue("displacement_absolute_tolerance",self.settings["displacement_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance",self.settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance",self.settings["residual_absolute_tolerance"])
        conv_params.AddValue("contact_displacement_relative_tolerance",self.contact_settings["contact_displacement_relative_tolerance"])
        conv_params.AddValue("contact_displacement_absolute_tolerance",self.contact_settings["contact_displacement_absolute_tolerance"])
        conv_params.AddValue("contact_residual_relative_tolerance",self.contact_settings["contact_residual_relative_tolerance"])
        conv_params.AddValue("contact_residual_absolute_tolerance",self.contact_settings["contact_residual_absolute_tolerance"])
        conv_params.AddValue("mortar_type",self.contact_settings["mortar_type"])
        conv_params.AddValue("condn_convergence_criterion",self.contact_settings["condn_convergence_criterion"])
        conv_params.AddValue("fancy_convergence_criterion",self.contact_settings["fancy_convergence_criterion"])
        conv_params.AddValue("print_convergence_criterion",self.contact_settings["print_convergence_criterion"])
        conv_params.AddValue("ensure_contact",self.contact_settings["ensure_contact"])
        conv_params.AddValue("gidio_debug",self.contact_settings["gidio_debug"])
        import contact_convergence_criteria_factory
        convergence_criterion = contact_convergence_criteria_factory.convergence_criterion(conv_params)
        return convergence_criterion.mechanical_convergence_criterion
    
    def _create_builder_and_solver(self):
        if  self.contact_settings["mortar_type"].GetString() != "":
            linear_solver = self.get_linear_solver()
            if self.settings["block_builder"].GetBool():
                if self.settings["multi_point_constraints_used"].GetBool():
                    raise Exception("MPCs not compatible with contact")
                else:
                    builder_and_solver = ContactStructuralMechanicsApplication.ContactResidualBasedBlockBuilderAndSolver(linear_solver)
            else:
                raise Exception("Contact not compatible with EliminationBuilderAndSolver")
        else:
            builder_and_solver = super()._create_builder_and_solver()
            
        return builder_and_solver
    
    def _create_mechanical_solver(self):
        if  self.contact_settings["mortar_type"].GetString() != "":
            if self.settings["analysis_type"].GetString() == "linear":
                mechanical_solver = self._create_linear_strategy()
            else:
                if(self.settings["line_search"].GetBool()):
                    mechanical_solver = self._create_contact_line_search_strategy()
                else:
                    mechanical_solver = self._create_contact_newton_raphson_strategy()
        else:
            mechanical_solver = super()._create_mechanical_solver()
                    
        return mechanical_solver
    
    def _create_contact_line_search_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        mechanical_convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
        newton_parameters = KratosMultiphysics.Parameters("""{}""")
        return ContactStructuralMechanicsApplication.LineSearchContactStrategy(computing_model_part, 
                                                                               mechanical_scheme, 
                                                                               linear_solver, 
                                                                               mechanical_convergence_criterion, 
                                                                               builder_and_solver, 
                                                                               self.settings["max_iteration"].GetInt(), 
                                                                               self.settings["compute_reactions"].GetBool(), 
                                                                               self.settings["reform_dofs_at_each_step"].GetBool(), 
                                                                               self.settings["move_mesh_flag"].GetBool(),
                                                                               newton_parameters
                                                                               )
    def _create_contact_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        mechanical_convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
        newton_parameters = KratosMultiphysics.Parameters("""{}""")
        newton_parameters.AddValue("adaptative_strategy",self.contact_settings["adaptative_strategy"])
        newton_parameters.AddValue("split_factor",self.contact_settings["split_factor"])
        newton_parameters.AddValue("max_number_splits",self.contact_settings["max_number_splits"])
        return ContactStructuralMechanicsApplication.ResidualBasedNewtonRaphsonContactStrategy(computing_model_part, 
                                                                                               mechanical_scheme, 
                                                                                               linear_solver, 
                                                                                               mechanical_convergence_criterion, 
                                                                                               builder_and_solver, 
                                                                                               self.settings["max_iteration"].GetInt(), 
                                                                                               self.settings["compute_reactions"].GetBool(), 
                                                                                               self.settings["reform_dofs_at_each_step"].GetBool(), 
                                                                                               self.settings["move_mesh_flag"].GetBool(),
                                                                                               newton_parameters,
                                                                                               self.processes_list
                                                                                               )
