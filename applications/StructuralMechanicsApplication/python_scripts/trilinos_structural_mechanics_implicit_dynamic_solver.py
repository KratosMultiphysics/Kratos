from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as SolidMechanicsApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# MPI
import KratosMultiphysics.mpi as mpi
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.MetisApplication as MetisApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import trilinos_structural_mechanics_solver

def CreateSolver(main_model_part, custom_settings):
    return TrilinosImplicitMechanicalSolver(main_model_part, custom_settings)

class TrilinosImplicitMechanicalSolver(trilinos_structural_mechanics_solver.TrilinosMechanicalSolver):

    ##constructor. the constructor shall only take care of storing the settings
    ##and the pointer to the main_model part. This is needed since at the point of constructing the
    ##model part is still not filled and the variables are not yet allocated
    ##
    ##real construction shall be delayed to the function "Initialize" which
    ##will be called once the model is already filled
    def __init__(self, main_model_part, custom_settings):

        #TODO: shall obtain the computing_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "trilinos_structural_mechanics_implicit_dynamic_solver",
            "echo_level": 0,
            "buffer_size": 2,
            "solution_type": "Dynamic",
            "time_integration_method": "Implicit",
            "scheme_type": "Newmark",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_step": false,
            "line_search": false,
            "implex": false,
            "compute_reactions": false,
            "compute_contact_forces": false,
            "block_builder": true,
            "clear_storage": false,
            "component_wise": false,
            "move_mesh_flag": true,
            "convergence_criterion": "Residual_criteria",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 10,
            "linear_solver_settings":{
                "solver_type" : "Klu",
                "scaling": false
            },
            "bodies_list": [],
            "problem_domain_sub_model_part_list": ["solid_model_part"],
            "processes_sub_model_part_list": [""]

        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        #construct the linear solver
        import trilinos_linear_solver_factory
        self.linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        print("Construction of Implicit Mechanical MPI Solver finished")

    def AddVariables(self):

        super(TrilinosImplicitMechanicalSolver, self).AddVariables()

        if self.settings["rotation_dofs"].GetBool():
            # Add specific variables for the problem (rotation dofs)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_TORQUE)

    def Initialize(self):

        print("::[Mechanical MPI Solver]:: -START-")
        # Construct the communicator
        self.EpetraCommunicator = TrilinosApplication.CreateCommunicator()

        # Get the solid computing model part
        self.computing_model_part = self.GetComputingModelPart()

        # Builder and solver creation
        builder_and_solver = self._GetBuilderAndSolver(self.settings["component_wise"].GetBool(),
                                                       self.settings["block_builder"].GetBool())

        # Solution scheme creation
        mechanical_scheme = self._GetSolutionScheme(self.settings["scheme_type"].GetString(),
                                                    self.settings["component_wise"].GetBool(),
                                                    self.settings["compute_contact_forces"].GetBool())

        # Get the convergence criterion
        mechanical_convergence_criterion = self._GetConvergenceCriterion()

        # Mechanical solver creation
        self._CreateMechanicalSolver(mechanical_scheme,
                                     mechanical_convergence_criterion,
                                     builder_and_solver,
                                     self.settings["max_iteration"].GetInt(),
                                     self.settings["compute_reactions"].GetBool(),
                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                     self.settings["move_mesh_flag"].GetBool(),
                                     self.settings["component_wise"].GetBool(),
                                     self.settings["line_search"].GetBool(),
                                     self.settings["implex"].GetBool())

        # Set echo_level
        self.mechanical_solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Check if everything is assigned correctly
        self.Check();

        print("::[Mechanical MPI Solver]:: -END- ")


    #### Decomposed Newton-Raphson resolution functions ####

    def SolverInitialize(self):
        self.mechanical_solver.Initialize()

    def SolverInitializeSolutionStep(self):
        self.mechanical_solver.InitializeSolutionStep()

    def SolverPredict(self):
        self.mechanical_solver.Predict()

    def SolverSolveSolutionStep(self):
        self.mechanical_solver.SolveSolutionStep()

    def SolverFinalizeSolutionStep(self):
        self.mechanical_solver.FinalizeSolutionStep()

    #### Specific internal functions ####

    def _GetSolutionScheme(self, scheme_type, component_wise, compute_contact_forces):

        if(scheme_type == "Newmark") or (scheme_type == "Bossak"):
            if (scheme_type == "Newmark"):
                mechanical_scheme = TrilinosApplication.TrilinosResidualBasedBossakDisplacementScheme(0.0)
            else:
                alpha = self.settings["damp_factor_m"].GetDouble()
                mechanical_scheme = TrilinosApplication.TrilinosResidualBasedBossakDisplacementScheme(alpha)

        return mechanical_scheme
