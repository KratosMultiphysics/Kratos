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
            "thermal_solve_frequency"        : 1,
            "voronoi_tesselation_frequency"  : 1000,
	        "porosity_update_frequency"      : 1000,
            "compute_motion"                 : true,
            "compute_direct_conduction"      : false,
            "compute_indirect_conduction"    : false,
            "compute_convection"             : false,
            "compute_radiation"              : false,
            "compute_friction_heat"          : false,
            "compute_adjusted_contact"       : false,
            "direct_conduction_model"        : "batchelor_obrien",
            "indirect_conduction_model"      : "surrounding_layer",
            "nusselt_correlation"            : "sphere_hanz_marshall",
            "radiation_model"                : "continuum_zhou",
            "adjusted_contact_model"         : "zhou",
            "voronoi_method"                 : "tesselation",
	        "porosity_method"                : "average_alpha_shape",
            "min_conduction_distance"        : 0.0000000275,
            "max_conduction_distance"        : 1.0,
            "fluid_layer_thickness"          : 0.4,
            "isothermal_core_radius"         : 0.5,
            "max_radiation_distance"         : 2.0,
            "global_porosity"                : 0.0,
            "alpha_shape_parameter"          : 1.2,
            "integral_tolerance"             : 0.000001,
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

        if (DEM_parameters["solver_settings"]["strategy"].GetString() != "thermal_sphere_strategy"):
            raise Exception('DEM', '"strategy" not available.')

        # General options
        self.compute_motion_option = self.thermal_settings["compute_motion"].GetBool()

        # Frequencies
        self.thermal_solve_frequency       = self.thermal_settings["thermal_solve_frequency"].GetInt()
        self.voronoi_tesselation_frequency = self.thermal_settings["voronoi_tesselation_frequency"].GetInt()
        self.porosity_update_frequency     = self.thermal_settings["porosity_update_frequency"].GetInt()

        # Set booleans for active heat transfer mechanisms
        self.compute_direct_conduction_option   = self.thermal_settings["compute_direct_conduction"].GetBool()
        self.compute_indirect_conduction_option = self.thermal_settings["compute_indirect_conduction"].GetBool()
        self.compute_convection_option          = self.thermal_settings["compute_convection"].GetBool()
        self.compute_radiation_option           = self.thermal_settings["compute_radiation"].GetBool()
        self.compute_friction_heat_option       = self.thermal_settings["compute_friction_heat"].GetBool()
        self.compute_adjusted_contact_option    = self.thermal_settings["compute_adjusted_contact"].GetBool()

        # Set models for heat transfer mechanisms
        self.direct_conduction_model   = self.thermal_settings["direct_conduction_model"].GetString()
        self.indirect_conduction_model = self.thermal_settings["indirect_conduction_model"].GetString()
        self.nusselt_correlation       = self.thermal_settings["nusselt_correlation"].GetString()
        self.radiation_model           = self.thermal_settings["radiation_model"].GetString()
        self.adjusted_contact_model    = self.thermal_settings["adjusted_contact_model"].GetString()
        self.voronoi_method            = self.thermal_settings["voronoi_method"].GetString()
        self.porosity_method           = self.thermal_settings["porosity_method"].GetString()

        # Check models
        if (self.direct_conduction_model != "batchelor_obrien" and
            self.direct_conduction_model != "thermal_pipe"     and
            self.direct_conduction_model != "collisional"):
            raise Exception('DEM', 'Direct thermal conduction model \'' + self.direct_conduction_model + '\' is not implemented.')

        if (self.indirect_conduction_model != "surrounding_layer" and
            self.indirect_conduction_model != "voronoi_a"         and
            self.indirect_conduction_model != "voronoi_b"         and
            self.indirect_conduction_model != "vargas_mccarthy"):
            raise Exception('DEM', 'Indirect thermal conduction model \'' + self.indirect_conduction_model + '\' is not implemented.')

        if (self.nusselt_correlation != "sphere_hanz_marshall" and
            self.nusselt_correlation != "sphere_whitaker"      and
            self.nusselt_correlation != "sphere_gunn"          and
            self.nusselt_correlation != "sphere_li_mason"):
            raise Exception('DEM', 'Nusselt number correlation \'' + self.nusselt_correlation + '\' is not implemented.')
        
        if (self.radiation_model != "continuum_zhou" and
            self.radiation_model != "continuum_krause"):
            raise Exception('DEM', 'Thermal radiation model \'' + self.radiation_model + '\' is not implemented.')

        if (self.adjusted_contact_model != "zhou" and
            self.adjusted_contact_model != "lu"):
            raise Exception('DEM', 'Adjusted contact model \'' + self.adjusted_contact_model + '\' is not implemented.')

        if (self.voronoi_method != "tesselation" and
            self.voronoi_method != "posority"):
            raise Exception('DEM', 'Voronoi method \'' + self.voronoi_method + '\' is not implemented.')
        
        if (self.porosity_method != "global"              and
            self.porosity_method != "average_convex_hull" and
            self.porosity_method != "average_alpha_shape"):
            raise Exception('DEM', 'Porosity method \'' + self.porosity_method + '\' is not implemented.')
        
        # Set model parameters
        self.min_conduction_distance = self.thermal_settings["min_conduction_distance"].GetDouble()
        self.max_conduction_distance = self.thermal_settings["max_conduction_distance"].GetDouble()
        self.fluid_layer_thickness   = self.thermal_settings["fluid_layer_thickness"].GetDouble()
        self.isothermal_core_radius  = self.thermal_settings["isothermal_core_radius"].GetDouble()
        self.max_radiation_distance  = self.thermal_settings["max_radiation_distance"].GetDouble()
        self.global_porosity         = self.thermal_settings["global_porosity"].GetDouble()
        self.alpha_parameter         = self.thermal_settings["alpha_shape_parameter"].GetDouble()
        self.integral_tolerance      = self.thermal_settings["integral_tolerance"].GetDouble()
        
        # Check / adjust parameters
        if (self.thermal_solve_frequency <= 0):
            self.thermal_solve_frequency = 1
        if (self.voronoi_tesselation_frequency < 0):
            self.voronoi_tesselation_frequency = 0
        if (self.porosity_update_frequency < 0):
            self.porosity_update_frequency = 0
        if (self.min_conduction_distance <= 0):
            raise Exception('DEM', '"min_conduction_distance" must be positive.')
        if (self.max_conduction_distance < 0):
            self.max_conduction_distance = 0
        if (self.fluid_layer_thickness < 0):
            self.fluid_layer_thickness = 0
        if (self.isothermal_core_radius < 0):
            self.isothermal_core_radius = 0
        if (self.isothermal_core_radius > 1):
            self.isothermal_core_radius = 1
        if (self.max_radiation_distance < 0 ):
            self.max_radiation_distance = 0
        if (self.global_porosity < 0 or self.global_porosity >= 1):
            raise Exception('DEM', '"global_porosity" must be between zero and one.')
        if (self.alpha_parameter < 0):
            raise Exception('DEM', '"alpha_shape_parameter" must be positive.')
        if (self.integral_tolerance <= 0):
            raise Exception('DEM', '"integral_tolerance" must be positive.')

        if (self.compute_indirect_conduction_option         and
           (self.indirect_conduction_model == "voronoi_a"   or
            self.indirect_conduction_model == "voronoi_b")  and
            self.voronoi_method            == "tesselation"):
            self.compute_voronoi = True
        else:
            self.compute_voronoi = False

        if   (self.compute_indirect_conduction_option         and
             (self.indirect_conduction_model == "voronoi_a"   or
              self.indirect_conduction_model == "voronoi_b")  and
              self.voronoi_method            == "posority"    and
              self.porosity_method           != "global"):
              self.compute_porosity = True
        elif (self.compute_convection_option                  and
             (self.nusselt_correlation == "sphere_gunn"       or
              self.nusselt_correlation == "sphere_li_mason")  and
              self.porosity_method     != "global"):
              self.compute_porosity = True
        elif (self.compute_radiation_option                   and
              self.radiation_model == "continuum_zhou"        and
              self.porosity_method != "global"):
              self.compute_porosity = True
        else:
              self.compute_porosity = False

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

        # Check / adjust fluid properties
        if (self.fluid_density              <= 0 or
            self.fluid_viscosity            <= 0 or
            self.fluid_thermal_conductivity <= 0 or
            self.fluid_heat_capacity        <= 0):
            raise Exception('DEM', '"global_fluid_properties" must contain positive values for material properties.')
        
        # Set booleans for graph writing
        self.PostGraphParticleTempMin   = DEM_parameters["PostGraphParticleTempMin"].GetBool()
        self.PostGraphParticleTempMax   = DEM_parameters["PostGraphParticleTempMax"].GetBool()
        self.PostGraphParticleTempAvg   = DEM_parameters["PostGraphParticleTempAvg"].GetBool()
        self.PostGraphParticleTempDev   = DEM_parameters["PostGraphParticleTempDev"].GetBool()
        self.PostGraphModelTempAvg      = DEM_parameters["PostGraphModelTempAvg"].GetBool()
        self.PostGraphFluxContributions = DEM_parameters["PostGraphHeatFluxContributions"].GetBool()

    def AddAdditionalVariables(self, model_part, DEM_parameters):
        # Add general additional variables (currently empty)
        BaseExplicitStrategy.AddAdditionalVariables(self, model_part, DEM_parameters)

        # Add thermal variables
        model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        model_part.AddNodalSolutionStepVariable(HEATFLUX)

    def SetAdditionalVariablesAndOptions(self):
        # Set general additional variables (currently empty)
        BaseExplicitStrategy.SetAdditionalVariablesAndOptions(self)

        # General options
        self.spheres_model_part.ProcessInfo.SetValue(THERMAL_FREQUENCY, self.thermal_solve_frequency)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, MOTION_OPTION, self.compute_motion_option)
        
        temperature_dependent_radius = False
        for properties in self.spheres_model_part.Properties:
            if ((properties.Has(THERMAL_EXPANSION_COEFFICIENT) and properties[THERMAL_EXPANSION_COEFFICIENT] != 0) or
                (properties.HasTable(TEMPERATURE,THERMAL_EXPANSION_COEFFICIENT))):
                temperature_dependent_radius = True
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, TEMPERATURE_DEPENDENT_RADIUS_OPTION, temperature_dependent_radius)

        # Booleans for active heat transfer mechanisms
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, DIRECT_CONDUCTION_OPTION,   self.compute_direct_conduction_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, INDIRECT_CONDUCTION_OPTION, self.compute_indirect_conduction_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, CONVECTION_OPTION,          self.compute_convection_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, RADIATION_OPTION,           self.compute_radiation_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, FRICTION_HEAT_OPTION,       self.compute_friction_heat_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, ADJUSTED_CONTACT_OPTION,    self.compute_adjusted_contact_option)

        # Models for heat transfer mechanisms
        self.spheres_model_part.ProcessInfo.SetValue(DIRECT_CONDUCTION_MODEL,   self.direct_conduction_model)
        self.spheres_model_part.ProcessInfo.SetValue(INDIRECT_CONDUCTION_MODEL, self.indirect_conduction_model)
        self.spheres_model_part.ProcessInfo.SetValue(CONVECTION_MODEL,          self.nusselt_correlation)
        self.spheres_model_part.ProcessInfo.SetValue(RADIATION_MODEL,           self.radiation_model)
        self.spheres_model_part.ProcessInfo.SetValue(ADJUSTED_CONTACT_MODEL,    self.adjusted_contact_model)
        self.spheres_model_part.ProcessInfo.SetValue(VORONOI_METHOD,            self.voronoi_method)
        self.spheres_model_part.ProcessInfo.SetValue(POSORITY_METHOD,           self.porosity_method)

        # Model parameters
        self.spheres_model_part.ProcessInfo.SetValue(MIN_CONDUCTION_DISTANCE,    self.min_conduction_distance)
        self.spheres_model_part.ProcessInfo.SetValue(MAX_CONDUCTION_DISTANCE,    self.max_conduction_distance)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_LAYER_THICKNESS,      self.fluid_layer_thickness)
        self.spheres_model_part.ProcessInfo.SetValue(ISOTHERMAL_CORE_RADIUS,     self.isothermal_core_radius)
        self.spheres_model_part.ProcessInfo.SetValue(MAX_RADIATION_DISTANCE,     self.max_radiation_distance)
        self.spheres_model_part.ProcessInfo.SetValue(AVERAGE_POROSITY,           self.global_porosity)
        self.spheres_model_part.ProcessInfo.SetValue(ALPHA_SHAPE_PARAMETER,      self.alpha_parameter)
        self.spheres_model_part.ProcessInfo.SetValue(INTEGRAL_TOLERANCE,         self.integral_tolerance)

        # Global properties for interstitial/surrounding fluid 
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_DENSITY,              self.fluid_density)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_VISCOSITY,            self.fluid_viscosity)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_THERMAL_CONDUCTIVITY, self.fluid_thermal_conductivity)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_HEAT_CAPACITY,        self.fluid_heat_capacity)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_TEMPERATURE,          self.fluid_temperature)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_VELOCITY,             self.fluid_velocity)

        # Booleans for writing graphs
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, GRAPH_PARTICLE_TEMP_MIN,                self.PostGraphParticleTempMin)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, GRAPH_PARTICLE_TEMP_MAX,                self.PostGraphParticleTempMax)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, GRAPH_PARTICLE_TEMP_AVG,                self.PostGraphParticleTempAvg)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, GRAPH_PARTICLE_TEMP_DEV,                self.PostGraphParticleTempDev)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, GRAPH_MODEL_TEMP_AVG,                   self.PostGraphModelTempAvg)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, GRAPH_PARTICLE_HEAT_FLUX_CONTRIBUTIONS, self.PostGraphFluxContributions)

    def Initialize(self):
        # Base class initializer
        BaseExplicitStrategy.Initialize(self)

        # Set thermal properties of provided in SubModelParts data
        (self.cplusplus_strategy).InitializeThermalDataInSubModelParts()
        (self.cplusplus_strategy).InitializeGraphOutput()

    def Predict(self):
        if (self.compute_motion_option):
            BaseExplicitStrategy.Predict(self)
    
    def InitializeSolutionStep(self):
        if (self.compute_motion_option):
            BaseExplicitStrategy.InitializeSolutionStep(self)
        else:
            (self.cplusplus_strategy).InitializeSolutionStep()

        # Perform tesselation-dependent tasks (triangulation or tetrahedralization)
        update_voronoi  = self.IsTimeToUpdateVoronoi()
        update_porosity = self.IsTimeToUpdatePorosity()
        if (update_voronoi or update_porosity):
            (self.cplusplus_strategy).TesselationTasks(update_voronoi,update_porosity)
    
    def IsTimeToUpdateVoronoi(self):
        if (self.compute_voronoi):
            step = self.spheres_model_part.ProcessInfo[TIME_STEPS]
            freq = self.voronoi_tesselation_frequency
            return step == 1 or (freq != 0 and step%freq == 0)
        else:
            return False
    
    def IsTimeToUpdatePorosity(self):
        if (self.compute_porosity):
            step = self.spheres_model_part.ProcessInfo[TIME_STEPS]
            freq = self.porosity_update_frequency
            return step == 1 or (freq != 0 and step%freq == 0)
        else:
            return False
    
    def SolveSolutionStep(self):
        if (self.compute_motion_option):
            (self.cplusplus_strategy).SolveSolutionStep()
        else:
            (self.cplusplus_strategy).SolveSolutionStepStatic()
        
        if (self.spheres_model_part.ProcessInfo[TEMPERATURE_DEPENDENT_RADIUS_OPTION]):
            (self.cplusplus_strategy).SetSearchRadiiOnAllParticles(self.spheres_model_part, self.search_increment, 1.0)
        
        return True

    def FinalizeSolutionStep(self):
        BaseExplicitStrategy.FinalizeSolutionStep(self)
        (self.cplusplus_strategy).WriteGraphOutput()
