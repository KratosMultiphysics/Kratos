# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

class ThermallyDrivenActiveFluidFormulation:
    """Helper class to define formulation parameters."""
    def __init__(self,settings):
        self.element_name = None
        self.element_integrates_in_time = False
        self.element_has_nodal_properties = False
        self.historical_nodal_variables_list = []
        self.non_historical_nodal_variables_list = []
        self.process_data = {}

        self.condition_name = "ThermallyDrivenActiveFluidCondition"

        if settings.Has("element_type"):
            formulation = settings["element_type"].GetString()
            if formulation == "divergence_constrained_mixed_dln":
                self._SetUpDivergenceConstrainedMixedDLN(settings)
            else:
                formulation_list = ["divergence_constrained_mixed_dln"]
                err_msg = f"Wrong \'element_type\' : \'{formulation}\' provided. Available options are:\n"
                for elem in formulation_list:
                    err_msg += f"\t- {elem}\n"
                # raise RuntimeError(err_msg) #TODO: Turn this into an error once the derived solvers handle this properly
                KratosMultiphysics.Logger.PrintWarning("ThermallyDrivenActiveFluidSolver", err_msg)
        else:
            print(settings)
            raise RuntimeError("Argument \'element_type\' not found in formulation settings.")

    def SetProcessInfo(self,model_part):
        for variable,value in self.process_data.items():
            model_part.ProcessInfo[variable] = value

    def _SetUpDivergenceConstrainedMixedDLN(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "divergence_constrained_mixed_dln",
            "gamma_stability_coefficient": 0.0,
            "rho_parameter": 0.0,
            "lambda_parameter": 0.0,
            "sigma_buoyancy_parameter": 0.0,
            "gravity_unitary_direction": [0.0, 0.0, 0.0],
            "time_integration_theta": 0.0
        }""")

        self.element_name = 'ThermallyDrivenActiveFluid'
        self.condition_name = "ThermallyDrivenActiveFluidCondition"

        settings.ValidateAndAssignDefaults(default_settings)

        # set the nodal material properties flag
        self.element_has_nodal_properties = True
        self.non_historical_nodal_variables_list = [KratosMultiphysics.VISCOSITY, KratosMultiphysics.DENSITY, KratosMultiphysics.CONDUCTIVITY]

        # set the process info variables values
        self.process_data[KratosMultiphysics.STABILIZATION_FACTOR] = settings["gamma_stability_coefficient"].GetDouble()
        self.process_data[KratosMultiphysics.Y1] = settings["rho_parameter"].GetDouble()
        self.process_data[KratosMultiphysics.Y2] = settings["lambda_parameter"].GetDouble()
        self.process_data[KratosMultiphysics.YF] = settings["sigma_buoyancy_parameter"].GetDouble()
        self.process_data[KratosMultiphysics.GRAVITY] = settings["gravity_unitary_direction"].GetVector()
        self.process_data[KratosMultiphysics.TIME_INTEGRATION_THETA] = settings["time_integration_theta"].GetDouble()

def CreateSolver(model, custom_settings):
    return ThermallyDrivenActiveFluidSolver(model, custom_settings)

class ThermallyDrivenActiveFluidSolver(FluidSolver):

    @classmethod
    def GetDefaultParameters(cls):

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "thermally_driven_active_fluid",
            "model_part_name": "FluidModelPart",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "reorder": false
            },
            "enforce_element_and_conditions_replacement": true,
            "material_import_settings": {
                "materials_filename": ""
            },
            "formulation": {
                "element_type": "divergence_constrained_mixed_dln",
                "gamma_stability_coefficient": 0.0,
                "rho_parameter": 0.0,
                "lambda_parameter": 0.0,
                "sigma_buoyancy_parameter": 0.0,
                "gravity_unitary_direction": [0.0, 0.0, 0.0],
                "time_integration_theta": 0.0
            },
            "maximum_iterations": 10,
            "echo_level": 0,
            "compute_reactions": false,
            "analysis_type": "non_linear",
            "reform_dofs_at_each_step": true,
            "assign_neighbour_elements_to_conditions": true,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_velocity_laplacian_tolerance": 1e-3,
            "absolute_velocity_laplacian_tolerance": 1e-5,
            "relative_temperature_tolerance": 1e-3,
            "absolute_temperature_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "relative_phi_tolerance": 1e-3,
            "absolute_phi_tolerance": 1e-5,
            "velocity_relaxation": 0.9,
            "pressure_relaxation": 0.9,
            "consider_periodic_conditions": false,
            "linear_solver_settings"        : {
                "solver_type" : "amgcl"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping"                : {
                "time_step"           : 0.0
            },
            "time_scheme":"steady",
            "move_mesh_flag": false
        }""")

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):

        super().__init__(model,custom_settings)

        # Set up the auxiliary class with the formulation settings
        self._SetFormulation()

        # Update the default buffer size according to the selected time scheme
        self._SetTimeSchemeBufferSize()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of ThermallyDrivenActiveFluidSolver finished.")

    def AddVariables(self):
        ## Add base class variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_LAPLACIAN)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.HEAT_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSUREAUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)

        # Adding variables required by the formulation (this includes the nodal material properties)
        for variable in self.historical_nodal_variables_list:
            self.main_model_part.AddNodalSolutionStepVariable(variable)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver variables added correctly.")

    def AddDofs(self):
        dofs_to_add = []
        dofs_to_add.append("VELOCITY_X")
        dofs_to_add.append("VELOCITY_Y")
        dofs_to_add.append("VELOCITY_Z")
        dofs_to_add.append("VELOCITY_LAPLACIAN_X")
        dofs_to_add.append("VELOCITY_LAPLACIAN_Y")
        dofs_to_add.append("VELOCITY_LAPLACIAN_Z")
        dofs_to_add.append("TEMPERATURE")
        dofs_to_add.append("PRESSURE")
        dofs_to_add.append("PRESSUREAUX")
        KratosMultiphysics.VariableUtils.AddDofsList(dofs_to_add, self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver DOFs added correctly.")

    def Initialize(self):
        # If the solver requires an instance of the stabilized formulation class, set the process info variables
        if hasattr(self, 'formulation'):
            self.formulation.SetProcessInfo(self.GetComputingModelPart())

        # Construct and initialize the solution strategy
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        solution_strategy.Initialize()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def InitializeSolutionStep(self):
        # Perform the solver InitializeSolutionStep
        self._GetSolutionStrategy().InitializeSolutionStep()

    def SolveSolutionStep(self):
        # Call the base fluid solver to solve current time step
        is_converged = super().SolveSolutionStep()

        # If the P2-P1 element is used, postprocess the pressure and the pressureaux variables in the quadratic nodes for the visualization
        # Note that this must be done in here (not in the FinalizeSolutionStep) in case the SolveSolutionStep
        # is called in a non-linear outer loop (e.g. from the FSI or the CHT solvers)
        if self.element_name == "ThermallyDrivenActiveFluid":
            KratosCFD.AutoFEFluidAuxiliaryUtilities.PostprocessP2P1ContinuousScalarVariable(self.GetComputingModelPart(), KratosMultiphysics.PRESSURE)
            KratosCFD.AutoFEFluidAuxiliaryUtilities.PostprocessP2P1ContinuousScalarVariable(self.GetComputingModelPart(), KratosMultiphysics.PRESSUREAUX)

        return is_converged

    def _SetFormulation(self):
        self.formulation = ThermallyDrivenActiveFluidFormulation(self.settings["formulation"])
        self.element_name = self.formulation.element_name
        self.condition_name = self.formulation.condition_name
        self.element_integrates_in_time = self.formulation.element_integrates_in_time
        self.element_has_nodal_properties = self.formulation.element_has_nodal_properties
        self.historical_nodal_variables_list = self.formulation.historical_nodal_variables_list
        self.non_historical_nodal_variables_list = self.formulation.non_historical_nodal_variables_list

    def _SetTimeSchemeBufferSize(self):
        scheme_type = self.settings["time_scheme"].GetString()
        if scheme_type == "steady":
            self.min_buffer_size = 3
        else:
            msg  = "Unknown time_scheme option found in project parameters:\n"
            msg += "\"" + scheme_type + "\"\n"
            msg += "Accepted values are \"dln\".\n"
            raise Exception(msg)
        
    def _ComputeDeltaTime(self):
        # User-defined delta time
        delta_time = self.settings["time_stepping"]["time_step"].GetDouble()
        return delta_time

    def _ComputeInitialDeltaTime(self):
        # User-defined delta time
        initial_delta_time = self.settings["time_stepping"]["time_step"].GetDouble()
        return initial_delta_time

    def _SetNodalProperties(self):
        set_viscosity = KratosMultiphysics.VISCOSITY in self.non_historical_nodal_variables_list
        set_density = KratosMultiphysics.DENSITY in self.non_historical_nodal_variables_list
        set_conductivity = KratosMultiphysics.CONDUCTIVITY in self.non_historical_nodal_variables_list

        # Get density, dynamic viscostity and conductivity from the properties of the first element
        for el in self.main_model_part.Elements:
            # Get DENSITY from properties
            if set_density:
                density_value = el.Properties.GetValue(KratosMultiphysics.DENSITY)
                if density_value <= 0.0:
                    raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(density_value,el.Properties.Id))
            # Get VISCOSITY from properties
            if set_viscosity:
                viscosity_value = el.Properties.GetValue(KratosMultiphysics.VISCOSITY)
                if viscosity_value <= 0.0:
                    raise Exception("VISCOSITY set to {0} in Properties {1}, positive number expected.".format(viscosity_value,el.Properties.Id))
            # Get VISCOSITY from properties
            if set_conductivity:
                conductivity_value = el.Properties.GetValue(KratosMultiphysics.CONDUCTIVITY)
            break
        else:
            raise Exception("No fluid elements found in the main model part.")

        # Transfer the obtained properties to the nodes
        if set_viscosity:
            self.main_model_part.ProcessInfo[KratosMultiphysics.VISCOSITY] = viscosity_value
        if set_density:
            self.main_model_part.ProcessInfo[KratosMultiphysics.DENSITY] = density_value
        if set_conductivity:
            self.main_model_part.ProcessInfo[KratosMultiphysics.CONDUCTIVITY] = conductivity_value