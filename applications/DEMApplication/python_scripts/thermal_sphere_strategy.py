from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import KratosMultiphysics.DEMApplication.sphere_strategy as SolverStrategy
BaseExplicitStrategy = SolverStrategy.ExplicitStrategy

class ExplicitStrategy(BaseExplicitStrategy):

    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures):
        # Initialize base class
        BaseExplicitStrategy.__init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures)

        # Get thermal settings and assign default values
        default_settings = Parameters("""
        {
            "compute_direct_conduction"      : false,
            "direct_conduction_model"        : "batchelor_obrien",
            "compute_indirect_conduction"    : false,
            "indirect_conduction_model"      : "surrounding_layer",
            "compute_convection"             : false,
            "nusselt_correlation"            : "sphere_hanz_marshall",
            "global_fluid_properties"        : {
                "fluid_density"              : 1.0,
                "fluid_viscosity"            : 1.0,
                "fluid_thermal_conductivity" : 1.0,
                "fluid_heat_capacity"        : 1.0,
                "fluid_temperature"          : 0.0,
                "fluid_velocity_X"           : 0.0,
                "fluid_velocity_Y"           : 0.0,
                "fluid_velocity_Z"           : 0.0
            }
        }""")

        self.thermal_settings = DEM_parameters["thermal_settings"]
        self.thermal_settings.ValidateAndAssignDefaults(default_settings)

        # Set booleans for active heat transfer mechanisms
        self.compute_direct_conduction_option   = DEM_parameters["compute_direct_conduction"].GetBool()
        self.compute_indirect_conduction_option = DEM_parameters["compute_indirect_conduction"].GetBool()
        self.compute_convection_option          = DEM_parameters["compute_convection"].GetBool()

        # Set models for heat transfer mechanisms
        self.direct_conduction_model   = DEM_parameters["direct_conduction_model"].GetString()
        self.indirect_conduction_model = DEM_parameters["indirect_conduction_model"].GetString()
        self.nusselt_correlation       = DEM_parameters["nusselt_correlation"].GetString()

        # Check input models
        if (self.direct_conduction_model != "batchelor_obrien" and
            self.direct_conduction_model != "thermal_pipe"     and
            self.direct_conduction_model != "collisional"):
            raise Exception('DEM', 'Direct thermal conduction model \'' + self.direct_conduction_model + '\' is not implemented.')

        if (self.indirect_conduction_model != "surrounding_layer"):
            raise Exception('DEM', 'Indirect thermal conduction model \'' + self.indirect_conduction_model + '\' is not implemented.')

        if (self.nusselt_correlation != "sphere_hanz_marshall" and
            self.nusselt_correlation != "sphere_whitaker"):
            raise Exception('DEM', 'Nusselt number correlation \'' + self.nusselt_correlation + '\' is not implemented.')

        # Set global properties of interstitial/surrounding fluid
        self.fluid_density              = DEM_parameters["fluid_density"].GetDouble()
        self.fluid_viscosity            = DEM_parameters["fluid_viscosity"].GetDouble()
        self.fluid_thermal_conductivity = DEM_parameters["fluid_thermal_conductivity"].GetDouble()
        self.fluid_heat_capacity        = DEM_parameters["fluid_heat_capacity"].GetDouble()
        self.fluid_temperature          = DEM_parameters["fluid_temperature"].GetDouble()
        self.fluid_velocity             = Vector(3)
        self.fluid_velocity[0]          = DEM_parameters["fluid_velocity_X"].GetDouble()
        self.fluid_velocity[1]          = DEM_parameters["fluid_velocity_Y"].GetDouble()
        self.fluid_velocity[2]          = DEM_parameters["fluid_velocity_Z"].GetDouble()

    def AddAdditionalVariables(self, model_part, DEM_parameters):
        # Add general additional variables (currently empty)
        BaseExplicitStrategy.AddAdditionalVariables(self, model_part, DEM_parameters)

        # Add thermal variables
        model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        model_part.AddNodalSolutionStepVariable(HEATFLUX)

    def SetAdditionalVariablesAndOptions(self):
        # Set general additional variables (currently empty)
        BaseExplicitStrategy.SetAdditionalVariablesAndOptions(self)

        # Booleans for active heat transfer mechanisms
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, DIRECT_CONDUCTION_OPTION,   self.compute_direct_conduction_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, INDIRECT_CONDUCTION_OPTION, self.compute_indirect_conduction_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, CONVECTION_OPTION,          self.compute_convection_option)

        # Models for heat transfer mechanisms
        self.spheres_model_part.ProcessInfo.SetValue(DIRECT_CONDUCTION_MODEL,   self.direct_conduction_model)
        self.spheres_model_part.ProcessInfo.SetValue(INDIRECT_CONDUCTION_MODEL, self.indirect_conduction_model)
        self.spheres_model_part.ProcessInfo.SetValue(CONVECTION_MODEL,          self.nusselt_correlation)

        # Global properties for interstitial/surrounding fluid 
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_DENSITY,              self.fluid_density)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_VISCOSITY,            self.fluid_viscosity)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_THERMAL_CONDUCTIVITY, self.fluid_thermal_conductivity)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_HEAT_CAPACITY,        self.fluid_heat_capacity)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_TEMPERATURE,          self.fluid_temperature)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_VELOCITY,             self.fluid_velocity)
