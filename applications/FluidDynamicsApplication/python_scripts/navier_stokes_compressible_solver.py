# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

## Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

from KratosMultiphysics.FluidDynamicsApplication import check_and_prepare_model_process_fluid

def CreateSolver(model, custom_settings):
    return NavierStokesCompressibleSolver(model, custom_settings)

class NavierStokesCompressibleSolver(FluidSolver):

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "compressible_solver_from_defaults",
            "model_part_name": "",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "two_element_test",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": "FluidMaterials.json"
            },
            "maximum_iterations": 10,
            "echo_level": 1,
            "time_order": 2,
            "time_scheme": "bdf2",
            "compute_reactions": false,
            "analysis_type": "non_linear",
            "reform_dofs_at_each_step" : true,
            "consider_periodic_conditions" : false,
            "relative_tolerance" : 1e-3,
            "absolute_tolerance" : 1e-5,
            "linear_solver_settings"       : {
                "solver_type"         : "amgcl",
                "max_iteration"       : 200,
                "tolerance"           : 1e-7,
                "provide_coordinates" : false,
                "smoother_type"       : "ilu0",
                "krylov_type"         : "gmres",
                "coarsening_type"     : "aggregation",
                "scaling"             : true,
                "verbosity"           : 0
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "assign_neighbour_elements_to_conditions": true,
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : true,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01
            },
            "periodic": "periodic",
            "move_mesh_flag": false
        }""")

        default_settings.AddMissingParameters(super(NavierStokesCompressibleSolver, cls).GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass=True # To be removed eventually
        # TODO: DO SOMETHING IN HERE TO REMOVE THE "time_order" FROM THE DEFAULT SETTINGS BUT KEEPING THE BACKWARDS COMPATIBILITY
        super(NavierStokesCompressibleSolver,self).__init__(model,custom_settings)

        self.min_buffer_size = 3
        self.element_name = "CompressibleNavierStokes"
        self.condition_name = "LineCondition"
        self.element_integrates_in_time = True
        self.element_has_nodal_properties = False

        ## Set the element replace settings
        #self._SetCompressibleElementReplaceSettings()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesCompressibleSolver finished.")

    def AddVariables(self):

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENTUM)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TOTAL_ENERGY)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONDUCTIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SPECIFIC_HEAT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.HEAT_CAPACITY_RATIO)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H) ## ?
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA) ## ?
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)  #for momentum
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.REACTION_DENSITY)  #for momentum
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.REACTION_ENERGY)  #for momentum

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE) ## ?
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.Y_WALL) ## ?
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.KINEMATIC_VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY)

        # Post-process
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.MACH)  #for momentum

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Monolithic compressible fluid solver variables added correctly")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MOMENTUM_X, KratosMultiphysics.REACTION_X, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MOMENTUM_Y, KratosMultiphysics.REACTION_Y, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MOMENTUM_Z, KratosMultiphysics.REACTION_Z, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DENSITY, KratosFluid.REACTION_DENSITY, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.TOTAL_ENERGY, KratosFluid.REACTION_ENERGY, self.main_model_part)

    def Initialize(self):
        # Construct and set the solution strategy
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        solution_strategy.Initialize()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def InitializeSolutionStep(self):
        (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)
        self._GetSolutionStrategy().InitializeSolutionStep()


    def Solve(self):
        (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)
        self._GetSolutionStrategy().Solve()

    # def PrepareModelPart(self):
    #     super(NavierStokesCompressibleSolver,self).PrepareModelPart()
    #     if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
    #         self._ExecuteAfterReading()

    # def _ExecuteAfterReading(self):
    #     ## Replace element and conditions
    #     KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

    #     ## Check that the input read has the shape we like
    #     prepare_model_part_settings = KratosMultiphysics.Parameters("{}")
    #     prepare_model_part_settings.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
    #     prepare_model_part_settings.AddValue("skin_parts",self.settings["skin_parts"])

    #     check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, prepare_model_part_settings).Execute()

    def _CreateScheme(self):
        domain_size = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        # Cases in which the element manages the time integration
        if self.element_integrates_in_time:
            # Rotation utility for compressible Navier-Stokes formulations
            # A custom rotation util is required as the nodal DOFs differs from the standard incompressible case
            rotation_utility = KratosFluid.CompressibleElementRotationUtility(
                domain_size,
                KratosMultiphysics.SLIP)
            # "Fake" scheme for those cases in where the element manages the time integration
            # It is required to perform the nodal update once the current time step is solved
            scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(rotation_utility)
            # In case the BDF2 scheme is used inside the element, the BDF time discretization utility is required to update the BDF coefficients
            if (self.settings["time_scheme"].GetString() == "bdf2"):
                time_order = 2
                self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
            else:
                err_msg = "Requested elemental time scheme \"" + self.settings["time_scheme"].GetString()+ "\" is not available.\n"
                err_msg += "Available options are: \"bdf2\""
                raise Exception(err_msg)
        # Cases in which a time scheme manages the time integration
        else:
            err_msg = "Custom scheme creation is not allowed. Compressible Navier-Stokes elements manage the time integration internally."
            raise Exception(err_msg)
        return scheme

    def _CreateConvergenceCriterion(self):
        convergence_criterion = KratosMultiphysics.ResidualCriteria(
            self.settings["relative_tolerance"].GetDouble(),
            self.settings["absolute_tolerance"].GetDouble())
        convergence_criterion.SetEchoLevel(self.settings["echo_level"].GetInt())
        return convergence_criterion