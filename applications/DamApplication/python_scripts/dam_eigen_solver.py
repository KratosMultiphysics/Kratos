from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication as KratosExternal
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.PoromechanicsApplication as KratosPoro

def CreateSolver(main_model_part, custom_settings):
    return DamEigenSolver(main_model_part, custom_settings)

class DamEigenSolver():

    def __init__(self, main_model_part, custom_settings):

        self.main_model_part = main_model_part

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "dam_eigen_solver",
            "echo_level": 0,
            "buffer_size": 1,
            "solution_type": "Dynamic",
            "analysis_type": "Linear",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "eigensolver_settings":{
                "solver_type": "FEAST",
                "print_feast_output"          : true,
                "perform_stochastic_estimate" : false,
                "solve_eigenvalue_problem"    : true,
                "compute_modal_contribution"  : false,
                "lambda_min"                  : 0.0,
                "lambda_max"                  : 500.0,
                "search_dimension"            : 4
            },
            "problem_domain_sub_model_part_list": ["solid_model_part"],
            "processes_sub_model_part_list": [""]
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.compute_modal_contribution = self.settings["eigensolver_settings"]["compute_modal_contribution"].GetBool()
        self.settings["eigensolver_settings"].RemoveValue("compute_modal_contribution")

        # eigensolver_settings are validated/assigned in the linear_solver
        print("Construction of Dam Eigensolver finished")


    def AddVariables(self):

        # Add displacements
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        # Add dynamic variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        # Add reactions for the displacements
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        print("::[Dam EigenSolver]:: Variables ADDED")


    def GetMinimumBufferSize(self):
        return 2;


    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)

        print("::[Dam EigenSolver]:: DOF's ADDED")


    def ImportModelPart(self):

        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):

            # Read ModelPart
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)

            # Create computing_model_part, set constitutive law and buffer size
            self._ExecuteAfterReading()

        else:
            raise Exception("Other input options are not yet implemented.")

        print ("Model reading finished")


    def Initialize(self):

        self.eigensolver_settings = self.settings["eigensolver_settings"]
        solver_type = self.eigensolver_settings["solver_type"].GetString()
        solution_type = self.settings["solution_type"].GetString()

        if solver_type == "FEAST":
            self.linear_solver = KratosExternal.FEASTSolver(self.eigensolver_settings)
        else:
            raise Exception("solver_type is not yet implemented.")

        if solution_type == "Dynamic":
            self.scheme = KratosStructural.EigensolverDynamicScheme()
        else:
            raise Exception("solution_type is not yet implemented.")

        self.builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosStructural.EigensolverStrategy(
            self.main_model_part,
            self.scheme,
            self.builder_and_solver,
            self.compute_modal_contribution)


    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.computing_model_part_name)


    def Solve(self):
        self.solver.Solve()


    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

        #### Specific internal functions ####

    def _ExecuteAfterReading(self):

        self.computing_model_part_name = "mechanical_computing_domain"

        # Auxiliary Kratos parameters object to be called by the CheckAndPepareModelProcess
        aux_params = KratosMultiphysics.Parameters("{}")
        aux_params.AddEmptyValue("computing_model_part_name").SetString(self.computing_model_part_name)

        # CheckAndPrepareModelProcess creates the solid_computational_model_part
        from KratosMultiphysics.PoromechanicsApplication import check_and_prepare_model_process_poro
        check_and_prepare_model_process_poro.CheckAndPrepareModelProcess(self.main_model_part, aux_params).Execute()

        # Constitutive law import
        from KratosMultiphysics.DamApplication import dam_constitutive_law_utility
        dam_constitutive_law_utility.SetConstitutiveLaw(self.main_model_part)

        self.main_model_part.SetBufferSize( self.settings["buffer_size"].GetInt() )
        minimum_buffer_size = self.GetMinimumBufferSize()
        if(minimum_buffer_size > self.main_model_part.GetBufferSize()):
            self.main_model_part.SetBufferSize( minimum_buffer_size )
