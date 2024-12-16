# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class files
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_monolithic_solver import StabilizedFormulation
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_monolithic_solver import NavierStokesMonolithicSolver

class LowMachFormulation(StabilizedFormulation):
    """Helper class to define stabilization-dependent parameters."""
    def __init__(self,settings):
        self.element_name = None
        self.element_integrates_in_time = False
        self.element_has_nodal_properties = False
        self.historical_nodal_variables_list = []
        self.non_historical_nodal_variables_list = []
        self.process_data = {}

        #TODO: Keep this until the MonolithicWallCondition is removed to ensure backwards compatibility in solvers with no defined condition_name
        self.condition_name = "MonolithicWallCondition"

        if settings.Has("element_type"):
            formulation = settings["element_type"].GetString()
            if formulation == "qsvms":
                self._SetUpQSVMS(settings)
            else:
                formulation_list = ["qsvms"]
                err_msg = f"Wrong \'element_type\' : \'{formulation}\' provided. Available options are:\n"
                for elem in formulation_list:
                    err_msg += f"\t- {elem}\n"
                # raise RuntimeError(err_msg) #TODO: Turn this into an error once the derived solvers handle this properly
                KratosMultiphysics.Logger.PrintWarning("NavierStokesLowMachSolver", err_msg)
        else:
            print(settings)
            raise RuntimeError("Argument \'element_type\' not found in formulation settings.")

    def _SetUpQSVMS(self,settings):

        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "qsvms",
            "dynamic_tau": 1.0
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "LowMachNavierStokes"
        self.condition_name = "LowMachNavierStokesWallCondition"
        self.element_integrates_in_time = True

        # set the nodal material properties flag
        self.element_has_nodal_properties = True
        self.historical_nodal_variables_list = []
        self.non_historical_nodal_variables_list = []

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()

def CreateSolver(model, custom_settings):
    return NavierStokesLowMachSolver(model, custom_settings)

class NavierStokesLowMachSolver(NavierStokesMonolithicSolver):

    @classmethod
    def GetDefaultParameters(cls):

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "low_mach",
            "model_part_name": "FluidModelPart",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "formulation": {
                "element_type": "qsvms"
            },
            "maximum_iterations": 10,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "compute_reactions": false,
            "analysis_type": "non_linear",
            "reform_dofs_at_each_step": true,
            "assign_neighbour_elements_to_conditions": true,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"        : {
                "solver_type" : "skyline_lu_factorization"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : false,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01,
                "time_step"           : 0.0
            },
            "time_scheme":"bdf2",
            "move_mesh_flag": false,
            "thermodynamic_pressure_settings" : {
                "flow_type" : "open",
                "value" : 101325.0
            }
        }""")

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesLowMachSolver finished.")

    def AddVariables(self):
        ## Add base class variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.HEAT_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        # Adding variables required by the formulation (this includes the nodal material properties)
        for variable in self.historical_nodal_variables_list:
            self.main_model_part.AddNodalSolutionStepVariable(variable)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver variables added correctly.")

    def AddDofs(self):
        dofs_and_reactions_to_add = []
        dofs_and_reactions_to_add.append(["VELOCITY_X", "REACTION_X"])
        dofs_and_reactions_to_add.append(["VELOCITY_Y", "REACTION_Y"])
        dofs_and_reactions_to_add.append(["VELOCITY_Z", "REACTION_Z"])
        dofs_and_reactions_to_add.append(["TEMPERATURE", "REACTION_FLUX"])
        dofs_and_reactions_to_add.append(["PRESSURE", "REACTION_WATER_PRESSURE"])
        KratosMultiphysics.VariableUtils.AddDofsList(dofs_and_reactions_to_add, self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver DOFs added correctly.")

    def _SetFormulation(self):
        self.formulation = LowMachFormulation(self.settings["formulation"])
        self.element_name = self.formulation.element_name
        self.condition_name = self.formulation.condition_name
        self.element_integrates_in_time = self.formulation.element_integrates_in_time
        self.element_has_nodal_properties = self.formulation.element_has_nodal_properties
        self.historical_nodal_variables_list = self.formulation.historical_nodal_variables_list
        self.non_historical_nodal_variables_list = self.formulation.non_historical_nodal_variables_list

    def _SetThermodynamicPressureSettings(self):
        th_pres_settings = self.settings["thermodynamic_pressure_settings"]
        flow_type = th_pres_settings["flow_type"].GetString()
        if flow_type == "open":
            self._GetComputingModelPart().ProcessInfo[KratosMultiphysics.PRESSURE] = th_pres_settings["value"].GetDouble()
        elif flow_type == "closed":
            raise Exception("Not implemented yet.")
        elif flow_type == "closed_with_inflow_outflow":
            raise Exception("Not implemented yet.")
        else:
            err_msg = f"Wrong value for 'flow_type'. Provided value is '{flow_type}'. Available options are:\n"
            err_msg += "\t- 'open'\n"
            err_msg += "\t- 'closed'\n"
            err_msg += "\t- 'closed_with_inflow_outflow'\n"

    #FIXME: We should fix the issue with the rotations
    def _CreateScheme(self):
        # "Fake" scheme for those cases in where the element manages the time integration
        # It is required to perform the nodal update once the current time step is solved
        # domain_size = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        # block_size = domain_size + 2
        # scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(domain_size, block_size)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

        # In case the BDF2 scheme is used inside the element, the BDF time discretization utility is required to update the BDF coefficients
        if (self.settings["time_scheme"].GetString() == "bdf2"):
            time_order = 2
            self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
        else:
            err_msg = "Requested elemental time scheme \"" + self.settings["time_scheme"].GetString()+ "\" is not available.\n"
            err_msg += "Available option is: \"bdf2\"."
            raise Exception(err_msg)

        return scheme
