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
            "integral_tolerance"             : 0.000001,
            "min_conduction_distance"        : 0.0000000275,
            "max_conduction_distance"        : 1.0,
            "fluid_layer_thickness"          : 0.4,
            "compute_convection"             : false,
            "nusselt_correlation"            : "sphere_hanz_marshall",
            "compute_radiation"              : false,
            "radiation_model"                : "continuum_zhou",
            "radiation_radius_factor"        : 3.0,
            "global_porosity"                : 0.0,
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
        self.compute_direct_conduction_option   = self.thermal_settings["compute_direct_conduction"].GetBool()
        self.compute_indirect_conduction_option = self.thermal_settings["compute_indirect_conduction"].GetBool()
        self.compute_convection_option          = self.thermal_settings["compute_convection"].GetBool()
        self.compute_radiation_option           = self.thermal_settings["compute_radiation"].GetBool()

        # Set models for heat transfer mechanisms
        self.direct_conduction_model   = self.thermal_settings["direct_conduction_model"].GetString()
        self.indirect_conduction_model = self.thermal_settings["indirect_conduction_model"].GetString()
        self.nusselt_correlation       = self.thermal_settings["nusselt_correlation"].GetString()
        self.radiation_model           = self.thermal_settings["radiation_model"].GetString()

        # Check models
        if (self.direct_conduction_model != "batchelor_obrien" and
            self.direct_conduction_model != "thermal_pipe"     and
            self.direct_conduction_model != "collisional"):
            raise Exception('DEM', 'Direct thermal conduction model \'' + self.direct_conduction_model + '\' is not implemented.')

        if (self.indirect_conduction_model != "surrounding_layer" and
            self.indirect_conduction_model != "voronoi_a"         and
            self.indirect_conduction_model != "voronoi_b"):
            raise Exception('DEM', 'Indirect thermal conduction model \'' + self.indirect_conduction_model + '\' is not implemented.')

        if (self.nusselt_correlation != "sphere_hanz_marshall" and
            self.nusselt_correlation != "sphere_whitaker"):
            raise Exception('DEM', 'Nusselt number correlation \'' + self.nusselt_correlation + '\' is not implemented.')
        
        if (self.radiation_model != "continuum_zhou" and
            self.radiation_model != "continuum_krause"):
            raise Exception('DEM', 'Thermal radiation model \'' + self.radiation_model + '\' is not implemented.')

        # Set model parameters
        self.integral_tolerance      = self.thermal_settings["integral_tolerance"].GetDouble()
        self.min_conduction_distance = self.thermal_settings["min_conduction_distance"].GetDouble()
        self.max_conduction_distance = self.thermal_settings["max_conduction_distance"].GetDouble()
        self.fluid_layer_thickness   = self.thermal_settings["fluid_layer_thickness"].GetDouble()
        self.radiation_radius_factor = self.thermal_settings["radiation_radius_factor"].GetDouble()
        self.global_porosity         = self.thermal_settings["global_porosity"].GetDouble()

        # Set global properties of interstitial/surrounding fluid
        self.fluid_props                = self.thermal_settings["global_fluid_properties"]
        self.fluid_density              = self.fluid_props["fluid_density"].GetDouble()
        self.fluid_viscosity            = self.fluid_props["fluid_viscosity"].GetDouble()
        self.fluid_thermal_conductivity = self.fluid_props["fluid_thermal_conductivity"].GetDouble()
        self.fluid_heat_capacity        = self.fluid_props["fluid_heat_capacity"].GetDouble()
        self.fluid_temperature          = self.fluid_props["fluid_temperature"].GetDouble()
        self.fluid_velocity             = Vector(3)
        self.fluid_velocity[0]          = self.fluid_props["fluid_velocity_X"].GetDouble()
        self.fluid_velocity[1]          = self.fluid_props["fluid_velocity_Y"].GetDouble()
        self.fluid_velocity[2]          = self.fluid_props["fluid_velocity_Z"].GetDouble()

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
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, RADIATION_OPTION,           self.compute_radiation_option)

        # Models for heat transfer mechanisms
        self.spheres_model_part.ProcessInfo.SetValue(DIRECT_CONDUCTION_MODEL,   self.direct_conduction_model)
        self.spheres_model_part.ProcessInfo.SetValue(INDIRECT_CONDUCTION_MODEL, self.indirect_conduction_model)
        self.spheres_model_part.ProcessInfo.SetValue(CONVECTION_MODEL,          self.nusselt_correlation)
        self.spheres_model_part.ProcessInfo.SetValue(RADIATION_MODEL,           self.radiation_model)

        # Model parameters
        self.spheres_model_part.ProcessInfo.SetValue(INTEGRAL_TOLERANCE,         self.integral_tolerance)
        self.spheres_model_part.ProcessInfo.SetValue(MIN_CONDUCTION_DISTANCE,    self.min_conduction_distance)
        self.spheres_model_part.ProcessInfo.SetValue(MAX_CONDUCTION_DISTANCE,    self.max_conduction_distance)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_LAYER_THICKNESS,      self.fluid_layer_thickness)
        self.spheres_model_part.ProcessInfo.SetValue(RADIATION_RADIUS_FACTOR,    self.radiation_radius_factor)
        self.spheres_model_part.ProcessInfo.SetValue(PRESCRIBED_GLOBAL_POROSITY, self.global_porosity)

        # Global properties for interstitial/surrounding fluid 
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_DENSITY,              self.fluid_density)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_VISCOSITY,            self.fluid_viscosity)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_THERMAL_CONDUCTIVITY, self.fluid_thermal_conductivity)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_HEAT_CAPACITY,        self.fluid_heat_capacity)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_TEMPERATURE,          self.fluid_temperature)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_VELOCITY,             self.fluid_velocity)

    def Initialize(self):
        # Base class initializer
        BaseExplicitStrategy.Initialize(self)

        # Set thermal properties of SubModelParts
        (self.cplusplus_strategy).InitializeThermalDataInSubModelParts()
