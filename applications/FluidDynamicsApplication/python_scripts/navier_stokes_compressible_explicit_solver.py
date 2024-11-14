# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

## Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
from KratosMultiphysics.FluidDynamicsApplication import check_and_prepare_model_process_fluid

def CreateSolver(model, custom_settings):
    return NavierStokesCompressibleExplicitSolver(model, custom_settings)

class NavierStokesCompressibleExplicitSolver(FluidSolver):
    def __init__(self, model, custom_settings):
        # Deprecated shock capturing
        self._CheckDeprecatedSettings(custom_settings)

        # Call base fluid solver constructor
        self._validate_settings_in_baseclass=True # To be removed eventually
        super(NavierStokesCompressibleExplicitSolver,self).__init__(model,custom_settings)

        # Define the formulation settings
        self.element_name = "CompressibleNavierStokesExplicit"
        if custom_settings["domain_size"].GetInt() == 2:
            self.condition_name = "LineCondition" # TODO: We need to create a Compressible NS condition (now using the base ones)
        elif custom_settings["domain_size"].GetInt() == 3:
            self.condition_name = "SurfaceCondition" # TODO: We need to create a Compressible NS condition (now using the base ones)
        else:
            err_msg = "Wrong domain size "
            raise Exception(err_msg)
        self.min_buffer_size = 2
        self.element_has_nodal_properties = False # Note that DENSITY is nodally stored but considered as a DOF

        KratosMultiphysics.Logger.PrintInfo("::[NavierStokesCompressibleExplicitSolver]:: ","Construction of NavierStokesCompressibleExplicitSolver finished.")

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "compressible_solver_from_defaults",
            "model_part_name": "FluidModelPart",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": "FluidMaterials.json"
            },
            "echo_level": 1,
            "time_order": 2,
            "time_scheme" : "RK4",
            "move_mesh_flag": false,
            "shock_capturing_settings" : { },
            "compute_reactions": false,
            "reform_dofs_at_each_step" : false,
            "assign_neighbour_elements_to_conditions": true,
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : true,
                "CFL_number"          : 1.0,
                "minimum_delta_time"  : 1.0e-8,
                "maximum_delta_time"  : 1.0e-2
            },
            "use_oss" : true
        }""")

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def AddVariables(self):
        # Add DOF variables (formulation written in conservative form) and reactions
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY) # Density DOF
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENTUM) # Momentum DOF
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TOTAL_ENERGY) # Total energy DOF
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.REACTION_DENSITY) # Density DOF reaction
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION) # Momentum DOF reaction
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.REACTION_ENERGY) # Total energy DOF reaction

        # Required variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.MASS_SOURCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.HEAT_SOURCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.NUMERICAL_ENTROPY) # TODO: This is only necessary whith shock capturing entropy_based
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.DENSITY_GRADIENT)

        # Post-process variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.MACH)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)

        KratosMultiphysics.Logger.PrintInfo("::[NavierStokesCompressibleExplicitSolver]:: ","Explicit compressible fluid solver variables added correctly")

    def AddDofs(self):
        dofs_with_reactions_list = []
        dofs_with_reactions_list.append(["DENSITY","REACTION_DENSITY"])
        dofs_with_reactions_list.append(["MOMENTUM_X","REACTION_X"])
        dofs_with_reactions_list.append(["MOMENTUM_Y","REACTION_Y"])
        if self.settings["domain_size"].GetInt() == 3:
            dofs_with_reactions_list.append(["MOMENTUM_Z","REACTION_Z"])
        dofs_with_reactions_list.append(["TOTAL_ENERGY","REACTION_ENERGY"])
        KratosMultiphysics.VariableUtils.AddDofsList(dofs_with_reactions_list, self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver DOFs added correctly.")

    def Initialize(self):
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.OSS_SWITCH] = int(self.settings["use_oss"].GetBool())
        self._ReadShockCapturingSettings()

        self.solver = self._get_solution_strategy()
        self.solver.SetEchoLevel(self.settings["echo_level"].GetInt())
        self.solver.Initialize()

        KratosMultiphysics.Logger.PrintInfo("::[NavierStokesCompressibleExplicitSolver]:: ","Explicit compressible fluid solver initialization finished.")

    def _get_solution_strategy(self):
        if not hasattr(self, '_solution_strategy'):
            self._solution_strategy = self._create_solution_strategy()
        return self._solution_strategy

    def _create_solution_strategy(self):
        self.computing_model_part = self.GetComputingModelPart()
        strategy_settings = KratosMultiphysics.Parameters('''{}''')
        strategy_settings.AddEmptyValue("rebuild_level").SetInt(1 if self.settings["reform_dofs_at_each_step"].GetBool() else 0)
        strategy_settings.AddEmptyValue("move_mesh_flag").SetBool(self.settings["move_mesh_flag"].GetBool())
        strategy_settings.AddEmptyValue("shock_capturing_settings").RecursivelyAddMissingParameters(self.settings["shock_capturing_settings"])

        requested_strategy = self.settings["time_scheme"].GetString()

        available_strategies = {
            "RK3-TVD"       : KratosFluid.CompressibleNavierStokesExplicitSolvingStrategyRungeKutta3TVD,
            "RK4"           : KratosFluid.CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4,
            "forward_euler" : KratosFluid.CompressibleNavierStokesExplicitSolvingStrategyForwardEuler,
            "bfecc" :         KratosFluid.CompressibleNavierStokesExplicitSolvingStrategyBFECC
        }

        if requested_strategy in available_strategies:
            strat = available_strategies[requested_strategy](self.computing_model_part, strategy_settings)
            self.settings["shock_capturing_settings"] = strategy_settings["shock_capturing_settings"]
            return strat

        err_msg = "Time scheme of type '{}' not available. Try any of\n".format(requested_strategy)
        for key in available_strategies:
            err_msg = err_msg + " - {}\n".format(key)
        raise RuntimeError(err_msg)

    @classmethod
    def _OverrideBoolParameterWithWarning(cls, parent, child, value):
        if parent.Has(child) and parent[child].GetBool() != value:
            KratosMultiphysics.Logger.PrintWarning("", "User-specifed {} will be overriden with {}".format(child, value))
        else:
            parent.AddEmptyValue(child)
        parent[child].SetBool(True)

    def _CreateEstimateDtUtility(self):
        """This method overloads FluidSolver in order to enforce:
        ```
        self.settings["time_stepping"]["consider_compressibility_in_CFL"] == True
        self.settings["time_stepping"]["nodal_density_formulation"] == True
        self.settings["time_stepping"]["consider_artificial_diffusion"] == SHOCK_CAPTURING_SWITCH
        ```
        """
        self._OverrideBoolParameterWithWarning(self.settings["time_stepping"], "consider_compressibility_in_CFL", True)
        self._OverrideBoolParameterWithWarning(self.settings["time_stepping"], "nodal_density_formulation", True)

        sc_enabled = self.GetComputingModelPart().ProcessInfo[KratosFluid.SHOCK_CAPTURING_SWITCH]
        self._OverrideBoolParameterWithWarning(self.settings["time_stepping"], "consider_artificial_diffusion", sc_enabled)

        estimate_dt_utility = KratosFluid.EstimateDtUtility(
                self.GetComputingModelPart(),
                self.settings["time_stepping"])

        return estimate_dt_utility

    def _CheckDeprecatedSettings(self, custom_settings):
        if custom_settings.Has("shock_capturing"):
            KratosMultiphysics.Logger.PrintInfo(
                '::[NavierStokesCompressibleExplicitSolver]:: "',
                "\n".join(['Setting "shock_capturing" is deprecated. Please use:',
                           '   {',
                           '       "shock_capturing_settings" : {',
                           '            "type" : "physics_based"',
                           '       }',
                           '   }',
                           'to maintain the same functionality with the newer syntax.'
                ]))

            custom_settings.RemoveValue("shock_capturing")
            # Not adding new syntax -> Using defauts


    def _ReadShockCapturingSettings(self):
        "Determines if shock capturing will be enabled and sets up SHOCK_CAPTURING_SWITCH"
        default_parameters = KratosMultiphysics.Parameters("""
            {
                    "type" : "physics_based",
                    "Parameters" : { }
            }
        """)

        self.settings["shock_capturing_settings"].ValidateAndAssignDefaults(default_parameters)
        enable_shock_capturing = self.settings["shock_capturing_settings"]["type"].GetString() != "none"

        self.GetComputingModelPart().ProcessInfo[KratosFluid.SHOCK_CAPTURING_SWITCH] = int(enable_shock_capturing)
