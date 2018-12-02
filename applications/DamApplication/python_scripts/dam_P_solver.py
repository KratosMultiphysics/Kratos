from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
# import KratosMultiphysics.ExternalSolversApplication as KratosExternal
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.DamApplication as KratosDam
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
#check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return DamUPSolver(main_model_part, custom_settings)

class DamUPSolver:
    def __init__(self, model_part, custom_settings):

        self.main_model_part = model_part

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type"                   : "dam_P_solver",
            "model_import_settings"         : {
                "input_type"       : "mdpa",
                "input_filename"   : "unknown_name",
                "input_file_label" : 0
            },
            "echo_level"                    : 1,
            "buffer_size"                   : 2,
            "processes_sub_model_part_list" : [""],
            "acoustic_solver_settings"      : {
                "strategy_type"               : "Newton-Raphson",
                "scheme_type"                 : "Newmark",
                "convergence_criterion"       : "Residual_criterion",
                "residual_relative_tolerance" : 0.0001,
                "residual_absolute_tolerance" : 1e-9,
                "max_iteration"               : 10,
                "move_mesh_flag"              : false,
                "echo_level"                  : 0,
                "linear_solver_settings"      : {
                    "solver_type"   : "AMGCL",
                    "max_iteration" : 200,
                    "tolerance"     : 1e-7,
                    "verbosity"     : 0,
                    "GMRES_size"    : 50
                    }
            }
         }""")

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        ## Definition of the linear solver
        amgcl_smoother = KratosMultiphysics.AMGCLSmoother.ILU0
        amgcl_krylov_type = KratosMultiphysics.AMGCLIterativeSolverType.BICGSTAB
        tolerance = self.settings["acoustic_solver_settings"]["linear_solver_settings"]["tolerance"].GetDouble()
        max_iterations = self.settings["acoustic_solver_settings"]["linear_solver_settings"]["max_iteration"].GetInt()
        verbosity = self.settings["acoustic_solver_settings"]["linear_solver_settings"]["verbosity"].GetInt()
        gmres_size = self.settings["acoustic_solver_settings"]["linear_solver_settings"]["GMRES_size"].GetInt()

        self.linear_solver =  KratosMultiphysics.AMGCLSolver(amgcl_smoother,amgcl_krylov_type,tolerance,max_iterations,verbosity,gmres_size)


    def AddVariables(self):
        ## Add pressure
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        ## Add Dynamic pressure Variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.Dt_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.Dt2_PRESSURE)

        print("Variables correctly added")


    def AddDofs(self):
        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.PRESSURE)

        print("P DOFs correctly added")

    def Initialize(self):

        if (self.settings["acoustic_solver_settings"]["scheme_type"].GetString() == "Newmark"):
            beta=0.25
            gamma=0.5
        else:
            raise Exception("Please use the Newmark Scheme, it is the only one available")

        move_mesh_flag = self.settings["acoustic_solver_settings"]["move_mesh_flag"].GetBool()
        max_iterations = self.settings["acoustic_solver_settings"]["max_iteration"].GetInt()
        res_rel_tol = self.settings["acoustic_solver_settings"]["residual_relative_tolerance"].GetDouble()
        res_abs_tol = self.settings["acoustic_solver_settings"]["residual_absolute_tolerance"].GetDouble()

        time_scheme = KratosDam.DamPScheme(beta,gamma)
        conv_criteria = KratosMultiphysics.ResidualCriteria(res_rel_tol,res_abs_tol)

        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(
            self.main_model_part,
            time_scheme,
            self.linear_solver,
            conv_criteria,
            max_iterations,
            False,
            False,
            move_mesh_flag)


        (self.solver).SetEchoLevel(self.settings["acoustic_solver_settings"]["echo_level"].GetInt())

        print ("Initialization Solver finished")

    def ImportModelPart(self):

        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):

            # Here it would be the place to import restart data if required
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)

            # Create computing_model_part, set constitutive law and buffer size
            self._ExecuteAfterReading()

        else:
            raise Exception("other input options are not yet implemented")

        print ("model reading finished")


    def GetMinimumBufferSize(self):
        return 2;

    def GetComputingModelPart(self):
        return self.main_model_part

    def Solve(self):
        (self.solver).Solve()

    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    def Clear(self):
        (self.solver).Clear()


    #### Specific internal functions ####

    def _ExecuteAfterReading(self):

        self.computing_model_part_name = "acoustic_computing_domain"

        # Auxiliary Kratos parameters object to be called by the CheckAndPepareModelProcess
        aux_params = KratosMultiphysics.Parameters("{}")
        aux_params.AddEmptyValue("computing_model_part_name").SetString(self.computing_model_part_name)

        # CheckAndPrepareModelProcess creates the solid_computational_model_part
        import check_and_prepare_model_process_poro
        check_and_prepare_model_process_poro.CheckAndPrepareModelProcess(self.main_model_part, aux_params).Execute()

        # Constitutive law import
        import dam_constitutive_law_utility
        dam_constitutive_law_utility.SetConstitutiveLaw(self.main_model_part)

        self.main_model_part.SetBufferSize( self.settings["buffer_size"].GetInt() )
        minimum_buffer_size = self.GetMinimumBufferSize()
        if(minimum_buffer_size > self.main_model_part.GetBufferSize()):
            self.main_model_part.SetBufferSize( minimum_buffer_size )

