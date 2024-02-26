# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.PoromechanicsApplication.poromechanics_U_Pl_solver import UPlSolver

def CreateSolver(model, custom_settings):
    return ExplicitUPlSolver(model, custom_settings)

class ExplicitUPlSolver(UPlSolver):
    """The Poromechanics explicit U (displacement) dynamic solver.

    This class creates the mechanical solvers for explicit dynamic analysis.
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super(ExplicitUPlSolver,self).__init__(model, custom_settings)

        self.min_buffer_size = 2           

        # Lumped mass-matrix is necessary for explicit analysis
        self.main_model_part.ProcessInfo[KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX] = True
        KratosMultiphysics.Logger.PrintInfo("::[ExplicitUPlSolver]:: Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "scheme_type"                : "Explicit_Central_Differences",
            "rebuild_level"              : 0,
            "theta_factor"               : 1.0,
            "g_factor"                   : 0.0,
            "calculate_xi"               : false,
            "xi_1_factor"                : 1.0
        }""")
        this_defaults.AddMissingParameters(super(ExplicitUPlSolver,self).GetDefaultParameters())
        return this_defaults

    def AddVariables(self):
        super(ExplicitUPlSolver,self).AddVariables()

        self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.DISPLACEMENT_OLD)
        # self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.DISPLACEMENT_OLDER)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.INTERNAL_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.FLUX_RESIDUAL)

        scheme_type = self.settings["scheme_type"].GetString()
        if(scheme_type == "Explicit_Velocity_Verlet"):
            self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.DAMPING_FORCE)
        
        # self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        # self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.RESIDUAL_VECTOR)

        KratosMultiphysics.Logger.PrintInfo("::[ExplicitUPlSolver]:: Variables ADDED")

    def AddDofs(self):
        # super(ExplicitUPlSolver,self).AddDofs()
        ## Solid dofs
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,self.main_model_part)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Z,self.main_model_part)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ACCELERATION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ACCELERATION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ACCELERATION_Z,self.main_model_part)

        KratosMultiphysics.VariableUtils().AddDof(KratosPoro.LIQUID_PRESSURE, KratosPoro.REACTION_LIQUID_PRESSURE,self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo("::[ExplicitUPlSolver]:: DOF's ADDED")

    def Initialize(self):
        # Using the base Initialize
        # super(ExplicitUPlSolver,self).Initialize()
        """Perform initialization after adding nodal variables and dofs to the main model part. """

        self.computing_model_part = self.GetComputingModelPart()

        # Fill the previous steps of the buffer with the initial conditions
        self._FillBuffer()

        # Solution scheme creation
        self.scheme = self._ConstructScheme(self.settings["scheme_type"].GetString())

        # Solver creation
        self.solver = self._ConstructSolver(self.settings["strategy_type"].GetString())

        # Set echo_level
        self.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Initialize Strategy
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        self.solver.Initialize()

        # Check if everything is assigned correctly
        self.Check()

        # Check and construct gp_to_nodal_variable process
        self._CheckAndConstructGPtoNodalVariableExtrapolationProcess()

    #### Specific internal functions ####
    def _ConstructScheme(self, scheme_type):
        scheme_type = self.settings["scheme_type"].GetString()

        # Setting the Rayleigh damping parameters
        process_info = self.main_model_part.ProcessInfo
        g_factor = self.settings["g_factor"].GetDouble()
        theta_factor = self.settings["theta_factor"].GetDouble()
        g_coeff = 0.0
        Dt = self.settings["time_step"].GetDouble()
        omega_1 = self.settings["omega_1"].GetDouble()
        omega_n = self.settings["omega_n"].GetDouble()
        rayleigh_alpha = self.settings["rayleigh_alpha"].GetDouble()
        rayleigh_beta = self.settings["rayleigh_beta"].GetDouble()
        if (scheme_type == "Explicit_Central_Differences" and g_factor >= 1.0):
            theta_factor = 0.5
            g_coeff = Dt*omega_n*omega_n*0.25*g_factor
        if self.settings["calculate_alpha_beta"].GetBool():
            xi_1 = self.settings["xi_1"].GetDouble()
            xi_n = self.settings["xi_n"].GetDouble()
            if (scheme_type == "Explicit_Central_Differences" and self.settings["calculate_xi"].GetBool()==True):
                xi_1_factor = self.settings["xi_1_factor"].GetDouble()                
                import numpy as np
                xi_1 = (np.sqrt(1+g_coeff*Dt)-theta_factor*omega_1*Dt*0.5)*xi_1_factor
                xi_n = (np.sqrt(1+g_coeff*Dt)-theta_factor*omega_n*Dt*0.5)
            rayleigh_beta = 2.0*(xi_n*omega_n-xi_1*omega_1)/(omega_n*omega_n-omega_1*omega_1)
            rayleigh_alpha = 2.0*xi_1*omega_1-rayleigh_beta*omega_1*omega_1
            KratosMultiphysics.Logger.PrintInfo("::[ExplicitUPlSolver]:: Scheme Information")
            KratosMultiphysics.Logger.PrintInfo("::[ExplicitUPlSolver]:: dt: ",Dt)
            KratosMultiphysics.Logger.PrintInfo("::[ExplicitUPlSolver]:: g_coeff: ",g_coeff)
            KratosMultiphysics.Logger.PrintInfo("::[ExplicitUPlSolver]:: omega_1: ",omega_1)
            KratosMultiphysics.Logger.PrintInfo("::[ExplicitUPlSolver]:: omega_n: ",omega_n)
            KratosMultiphysics.Logger.PrintInfo("::[ExplicitUPlSolver]:: xi_1: ",xi_1)
            KratosMultiphysics.Logger.PrintInfo("::[ExplicitUPlSolver]:: xi_n: ",xi_n)
            KratosMultiphysics.Logger.PrintInfo("::[ExplicitUPlSolver]:: Alpha and Beta output")
            KratosMultiphysics.Logger.PrintInfo("::[ExplicitUPlSolver]:: rayleigh_alpha: ",rayleigh_alpha)
            KratosMultiphysics.Logger.PrintInfo("::[ExplicitUPlSolver]:: rayleigh_beta: ",rayleigh_beta)
        
        process_info.SetValue(StructuralMechanicsApplication.RAYLEIGH_ALPHA, rayleigh_alpha)
        process_info.SetValue(StructuralMechanicsApplication.RAYLEIGH_BETA, rayleigh_beta)
        process_info.SetValue(KratosPoro.G_COEFFICIENT, g_coeff)
        process_info.SetValue(KratosPoro.THETA_FACTOR, theta_factor)

        # Setting the time integration schemes
        if(scheme_type == "Explicit_Central_Differences"):
            scheme = KratosPoro.PoroExplicitCDScheme()
        elif(scheme_type == "Explicit_Velocity_Verlet"):
            scheme = KratosPoro.PoroExplicitVVScheme()
        else:
            err_msg =  "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"Explicit_Central_Differences\", \"Explicit_Velocity_Verlet\" "
            raise Exception(err_msg)
        return scheme

    def _ConstructSolver(self, strategy_type):
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.ERROR_RATIO, self.settings["displacement_relative_tolerance"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.ERROR_INTEGRATION_POINT, self.settings["displacement_absolute_tolerance"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosPoro.IS_CONVERGED, True)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.NL_ITERATION_NUMBER, 1)

        nonlocal_damage = self.settings["nonlocal_damage"].GetBool()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        reform_step_dofs = self.settings["reform_dofs_at_each_step"].GetBool()
        move_mesh_flag = self.settings["move_mesh_flag"].GetBool()

        self.strategy_params = KratosMultiphysics.Parameters("{}")
        self.strategy_params.AddValue("loads_sub_model_part_list",self.loads_sub_sub_model_part_list)
        self.strategy_params.AddValue("loads_variable_list",self.settings["loads_variable_list"])
        # NOTE: A rebuild level of 0 means that the nodal mass is calculated only once at the beginning (Initialize)
        #       A rebuild level higher than 0 means that the nodal mass can be updated at the beginning of each step (InitializeSolutionStep)
        self.strategy_params.AddValue("rebuild_level",self.settings["rebuild_level"])

        if nonlocal_damage:
            self.strategy_params.AddValue("body_domain_sub_model_part_list",self.body_domain_sub_sub_model_part_list)
            self.strategy_params.AddValue("characteristic_length",self.settings["characteristic_length"])
            self.strategy_params.AddValue("search_neighbours_step",self.settings["search_neighbours_step"])
            solving_strategy = KratosPoro.PoromechanicsExplicitNonlocalStrategy(self.computing_model_part,
                                                                            self.scheme,
                                                                            self.strategy_params,
                                                                            compute_reactions,
                                                                            reform_step_dofs,
                                                                            move_mesh_flag)
        else:
            solving_strategy = KratosPoro.PoromechanicsExplicitStrategy(self.computing_model_part,
                                                                            self.scheme,
                                                                            self.strategy_params,
                                                                            compute_reactions,
                                                                            reform_step_dofs,
                                                                            move_mesh_flag)

        return solving_strategy

    #### Private functions ####

