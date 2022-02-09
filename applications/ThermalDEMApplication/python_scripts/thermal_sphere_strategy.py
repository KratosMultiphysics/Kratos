# Imports
from   KratosMultiphysics import *
from   KratosMultiphysics.DEMApplication import *
from   KratosMultiphysics.ThermalDEMApplication import *
import KratosMultiphysics.DEMApplication.sphere_strategy as SolverStrategy
import KratosMultiphysics.ThermalDEMApplication.default_input_settings as DefaultSettings

# Set base class
BaseStrategy = SolverStrategy.ExplicitStrategy

# Auxiliary functions
def GetBoolParameterIfItExists(parameters, key):
    if key in parameters.keys():
        return parameters[key].GetBool()
    else:
        return False

# Strategy class
class ExplicitStrategy(BaseStrategy):

    ####################################### DERIVED METHODS #######################################
    #----------------------------------------------------------------------------------------------
    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures):
        # Initialize base class
        BaseStrategy.__init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures)

        # Get and validate input parameters
        self.GetProjectParameters(DEM_parameters)
        self.CheckProjectParameters()

        # Set flags
        self.SetVoronoiPorosityFlags()
        self.SetGraphFlags()

        # Create utilities
        self.CreateCPlusPlusUtilities()

    #----------------------------------------------------------------------------------------------
    def CreateCPlusPlusStrategy(self):
        # Set standard options
        BaseStrategy.SetVariablesAndOptions(self)

        # Set thermal options
        self.SetThermalVariablesAndOptions()

        # Create cpp strategy object
        translational_integration_scheme = self.DEM_parameters["TranslationalIntegrationScheme"].GetString()
        
        if (translational_integration_scheme == 'Velocity_Verlet'):
            raise Exception('ThermalDEM', '"Thermal strategy for translational integration scheme \'' + translational_integration_scheme + '\' is not implemented.')
        else:
            self.cplusplus_strategy = ThermalExplicitSolverStrategy(self.settings,
                                                                    self.max_delta_time,
                                                                    self.n_step_search,
                                                                    self.safety_factor,
                                                                    self.delta_option,
                                                                    self.creator_destructor,
                                                                    self.dem_fem_search,
                                                                    self.search_strategy,
                                                                    self.solver_settings)

    #----------------------------------------------------------------------------------------------
    def ModifyProperties(self, properties, param = 0):
        if param:
            return

        # Set standard properties
        BaseStrategy.ModifyProperties(self, properties, param)

        # Set thermal integration scheme
        self.SetThermalIntegrationScheme(properties)
    
    #----------------------------------------------------------------------------------------------
    def AddVariables(self):
        # Add standard variables
        BaseStrategy.AddVariables(self)

        # Add thermal variables to all model parts
        self.spheres_model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        self.cluster_model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        self.inlet_model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        self.fem_model_part.AddNodalSolutionStepVariable(TEMPERATURE)

        self.spheres_model_part.AddNodalSolutionStepVariable(HEATFLUX)
        self.cluster_model_part.AddNodalSolutionStepVariable(HEATFLUX)
        self.inlet_model_part.AddNodalSolutionStepVariable(HEATFLUX)
        self.fem_model_part.AddNodalSolutionStepVariable(HEATFLUX)
    
    #----------------------------------------------------------------------------------------------
    def Initialize(self):
        # Base class initializer
        BaseStrategy.Initialize(self)

        # Initialize utilities
        self.InitializeCPlusPlusUtilities()
    
    #----------------------------------------------------------------------------------------------
    def InitializeSolutionStep(self):
        if (self.compute_motion_option):
            BaseStrategy.InitializeSolutionStep(self)
        else:
            (self.cplusplus_strategy).InitializeSolutionStep()

        # Perform tesselation-dependent tasks (triangulation or tetrahedralization)
        if (self.IsTimeToUpdateVoronoi() or self.IsTimeToUpdatePorosity()):
            self.tesselation_utils.ExecuteInitializeSolutionStep(self.spheres_model_part)
    
    #----------------------------------------------------------------------------------------------
    def Predict(self):
        if (self.compute_motion_option):
            BaseStrategy.Predict(self)
    
    #----------------------------------------------------------------------------------------------
    def SolveSolutionStep(self):
        # Solve step according to motion type
        if (self.compute_motion_option):
            (self.cplusplus_strategy).SolveSolutionStep()
        else:
            (self.cplusplus_strategy).SolveSolutionStepStatic()
        
        # Update temperature dependent search radii
        if (self.spheres_model_part.ProcessInfo[TEMPERATURE_DEPENDENT_RADIUS_OPTION]):
            (self.cplusplus_strategy).SetSearchRadiiOnAllParticles(self.spheres_model_part, self.search_increment, 1.0)
        
        return True
    
    #----------------------------------------------------------------------------------------------
    def FinalizeSolutionStep(self):
        BaseStrategy.FinalizeSolutionStep(self)

        # Write output graphs
        if (self.write_graph):
            self.graph_utils.ExecuteFinalizeSolutionStep(self.spheres_model_part)

    #----------------------------------------------------------------------------------------------
    def Finalize(self):
        BaseStrategy.Finalize(self)

        # Close graph files
        if (self.write_graph):
            self.graph_utils.ExecuteFinalize()

    ####################################### PARTICULAR METHODS #######################################
    #----------------------------------------------------------------------------------------------
    def GetProjectParameters(self, DEM_parameters):
        # Get thermal settings and assign default values (in case it was not previously done)
        default_settings = DefaultSettings.GetDefaultInputSettings()

        if "thermal_settings" in self.DEM_parameters.keys():
            self.thermal_settings = DEM_parameters["thermal_settings"]
        else:
            self.thermal_settings = Parameters("""{}""")
        
        self.thermal_settings.ValidateAndAssignDefaults(default_settings)

        # General options
        self.compute_motion_option = self.thermal_settings["compute_motion"].GetBool()

        # Integration scheme
        self.thermal_integration_scheme = self.thermal_settings["ThermalIntegrationScheme"].GetString()

        # Frequencies
        self.thermal_solve_frequency       = self.thermal_settings["thermal_solve_frequency"].GetInt()
        self.voronoi_tesselation_frequency = self.thermal_settings["voronoi_tesselation_frequency"].GetInt()
        self.porosity_update_frequency     = self.thermal_settings["porosity_update_frequency"].GetInt()

        # Active heat transfer mechanisms
        self.compute_direct_conduction_option   = GetBoolParameterIfItExists(self.thermal_settings, "compute_direct_conduction")
        self.compute_indirect_conduction_option = GetBoolParameterIfItExists(self.thermal_settings, "compute_indirect_conduction")
        self.compute_convection_option          = GetBoolParameterIfItExists(self.thermal_settings, "compute_convection")
        self.compute_radiation_option           = GetBoolParameterIfItExists(self.thermal_settings, "compute_radiation")
        self.compute_friction_heat_option       = GetBoolParameterIfItExists(self.thermal_settings, "compute_friction_heat")
        self.compute_adjusted_contact_option    = GetBoolParameterIfItExists(self.thermal_settings, "compute_adjusted_contact")

        # Models for heat transfer
        self.direct_conduction_model   = self.thermal_settings["direct_conduction_model"].GetString()
        self.indirect_conduction_model = self.thermal_settings["indirect_conduction_model"].GetString()
        self.nusselt_correlation       = self.thermal_settings["nusselt_correlation"].GetString()
        self.radiation_model           = self.thermal_settings["radiation_model"].GetString()
        self.adjusted_contact_model    = self.thermal_settings["adjusted_contact_model"].GetString()
        self.voronoi_method            = self.thermal_settings["voronoi_method"].GetString()
        self.porosity_method           = self.thermal_settings["porosity_method"].GetString()
        
        # Model parameters
        self.min_conduction_distance  = self.thermal_settings["min_conduction_distance"].GetDouble()
        self.max_conduction_distance  = self.thermal_settings["max_conduction_distance"].GetDouble()
        self.fluid_layer_thickness    = self.thermal_settings["fluid_layer_thickness"].GetDouble()
        self.isothermal_core_radius   = self.thermal_settings["isothermal_core_radius"].GetDouble()
        self.max_radiation_distance   = self.thermal_settings["max_radiation_distance"].GetDouble()
        self.friction_heat_conversion = self.thermal_settings["friction_heat_conversion_ratio"].GetDouble()
        self.global_porosity          = self.thermal_settings["global_porosity"].GetDouble()
        self.alpha_parameter          = self.thermal_settings["alpha_shape_parameter"].GetDouble()
        self.integral_tolerance       = self.thermal_settings["integral_tolerance"].GetDouble()

        # Interstitial fluid properties
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
        
        # Graph writing
        self.PostGraphParticleTempMin   = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphParticleTempMin")
        self.PostGraphParticleTempMax   = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphParticleTempMax")
        self.PostGraphParticleTempAvg   = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphParticleTempAvg")
        self.PostGraphParticleTempDev   = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphParticleTempDev")
        self.PostGraphModelTempAvg      = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphModelTempAvg")
        self.PostGraphFluxContributions = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphFluxContributions")

    #----------------------------------------------------------------------------------------------
    def CheckProjectParameters(self):
        # Integration scheme
        if self.thermal_integration_scheme == "Forward_Euler":
            self.thermal_integration_scheme = "ThermalForwardEulerScheme"
        else:
            raise Exception('ThermalDEM', 'Integration scheme \'' + self.thermal_integration_scheme + '\' is not implemented.')
        
        # Models for heat transfer
        if (self.direct_conduction_model != "batchelor_obrien" and
            self.direct_conduction_model != "thermal_pipe"     and
            self.direct_conduction_model != "collisional"):
            raise Exception('ThermalDEM', 'Direct thermal conduction model \'' + self.direct_conduction_model + '\' is not implemented.')

        if (self.indirect_conduction_model != "surrounding_layer" and
            self.indirect_conduction_model != "voronoi_a"         and
            self.indirect_conduction_model != "voronoi_b"         and
            self.indirect_conduction_model != "vargas_mccarthy"):
            raise Exception('ThermalDEM', 'Indirect thermal conduction model \'' + self.indirect_conduction_model + '\' is not implemented.')

        if (self.nusselt_correlation != "sphere_hanz_marshall" and
            self.nusselt_correlation != "sphere_whitaker"      and
            self.nusselt_correlation != "sphere_gunn"          and
            self.nusselt_correlation != "sphere_li_mason"):
            raise Exception('ThermalDEM', 'Nusselt number correlation \'' + self.nusselt_correlation + '\' is not implemented.')
        
        if (self.radiation_model != "continuum_zhou" and
            self.radiation_model != "continuum_krause"):
            raise Exception('ThermalDEM', 'Thermal radiation model \'' + self.radiation_model + '\' is not implemented.')

        if (self.adjusted_contact_model != "zhou" and
            self.adjusted_contact_model != "lu"   and
            self.adjusted_contact_model != "morris"):
            raise Exception('ThermalDEM', 'Adjusted contact model \'' + self.adjusted_contact_model + '\' is not implemented.')

        if (self.voronoi_method != "tesselation" and
            self.voronoi_method != "posority"):
            raise Exception('ThermalDEM', 'Voronoi method \'' + self.voronoi_method + '\' is not implemented.')
        
        if (self.porosity_method != "global"              and
            self.porosity_method != "average_convex_hull" and
            self.porosity_method != "average_alpha_shape"):
            raise Exception('ThermalDEM', 'Porosity method \'' + self.porosity_method + '\' is not implemented.')

        # Model parameters
        if (self.thermal_solve_frequency <= 0):
            self.thermal_solve_frequency = 1
        if (self.voronoi_tesselation_frequency < 0):
            self.voronoi_tesselation_frequency = 0
        if (self.porosity_update_frequency < 0):
            self.porosity_update_frequency = 0
        if (self.min_conduction_distance <= 0):
            raise Exception('ThermalDEM', '"min_conduction_distance" must be positive.')
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
        if (self.friction_heat_conversion < 0 or self.friction_heat_conversion > 1):
            raise Exception('ThermalDEM', '"friction_heat_conversion_ratio" must be between zero and one.')
        if (self.global_porosity < 0 or self.global_porosity >= 1):
            raise Exception('ThermalDEM', '"global_porosity" must be between zero and one.')
        if (self.alpha_parameter < 0):
            raise Exception('ThermalDEM', '"alpha_shape_parameter" must be positive.')
        if (self.integral_tolerance <= 0):
            raise Exception('ThermalDEM', '"integral_tolerance" must be positive.')
        
        # Fluid properties
        if (self.fluid_density              <= 0 or
            self.fluid_viscosity            <= 0 or
            self.fluid_thermal_conductivity <= 0 or
            self.fluid_heat_capacity        <= 0):
            raise Exception('ThermalDEM', '"global_fluid_properties" must contain positive values for material properties.')

    #----------------------------------------------------------------------------------------------
    def SetVoronoiPorosityFlags(self):
        # Flag for computing voronoi diagram in a given frequency
        if (self.compute_indirect_conduction_option         and
           (self.indirect_conduction_model == "voronoi_a"   or
            self.indirect_conduction_model == "voronoi_b")  and
            self.voronoi_method            == "tesselation"):
            self.compute_voronoi = True
        else:
            self.compute_voronoi = False

        # Flag for computing porosity in a given frequency
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
    
    #----------------------------------------------------------------------------------------------
    def SetGraphFlags(self):
        if (self.PostGraphParticleTempMin  or
            self.PostGraphParticleTempMax  or
            self.PostGraphParticleTempAvg  or 
            self.PostGraphParticleTempDev  or
            self.PostGraphModelTempAvg     or
            self.PostGraphFluxContributions):
            self.write_graph = True
        else:
            self.write_graph = False

    #----------------------------------------------------------------------------------------------
    def CreateCPlusPlusUtilities(self):
        self.thermal_data_utils = SetThermalDataUtilities()

        if (self.compute_voronoi or self.compute_porosity):
            if self.dimension == 2:
                self.tesselation_utils = TesselationUtilities2D()
            elif self.dimension == 3:
                self.tesselation_utils = TesselationUtilities3D()
            
        if (self.write_graph):
            self.graph_utils = GraphUtilities()
    
    #----------------------------------------------------------------------------------------------
    def SetThermalIntegrationScheme(self, properties):
        integration_scheme_name = self.thermal_integration_scheme

        if properties.Has(DEM_THERMAL_INTEGRATION_SCHEME_NAME):
            if properties[DEM_THERMAL_INTEGRATION_SCHEME_NAME] == "Forward_Euler":
                integration_scheme_name = "ThermalForwardEulerScheme"
            else:
                raise Exception('ThermalDEM', 'Integration scheme \'' + integration_scheme_name + '\' is not implemented.')
        
        try:
            thermal_scheme = eval(integration_scheme_name)()
        except:
            raise Exception('The class corresponding to the thermal integration scheme named ' + integration_scheme_name + ' has not been added to python. Please, select a different name or add the required class.')
        
        thermal_scheme.SetThermalIntegrationSchemeInProperties(properties, True)

    #----------------------------------------------------------------------------------------------
    def SetThermalVariablesAndOptions(self):
        # General options
        self.spheres_model_part.ProcessInfo.SetValue(THERMAL_FREQUENCY, self.thermal_solve_frequency)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, MOTION_OPTION, self.compute_motion_option)
        
        temperature_dependent_radius = False
        for properties in self.spheres_model_part.Properties:
            if ((properties.Has(THERMAL_EXPANSION_COEFFICIENT) and properties[THERMAL_EXPANSION_COEFFICIENT] != 0) or
                (properties.HasTable(TEMPERATURE,THERMAL_EXPANSION_COEFFICIENT))):
                temperature_dependent_radius = True
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, TEMPERATURE_DEPENDENT_RADIUS_OPTION, temperature_dependent_radius)

        # Active heat transfer mechanisms
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, DIRECT_CONDUCTION_OPTION,   self.compute_direct_conduction_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, INDIRECT_CONDUCTION_OPTION, self.compute_indirect_conduction_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, CONVECTION_OPTION,          self.compute_convection_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, RADIATION_OPTION,           self.compute_radiation_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, FRICTION_HEAT_OPTION,       self.compute_friction_heat_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, ADJUSTED_CONTACT_OPTION,    self.compute_adjusted_contact_option)

        # Models for heat transfer
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
        self.spheres_model_part.ProcessInfo.SetValue(FRICTION_HEAT_CONVERSION,   self.friction_heat_conversion)
        self.spheres_model_part.ProcessInfo.SetValue(AVERAGE_POROSITY,           self.global_porosity)
        self.spheres_model_part.ProcessInfo.SetValue(ALPHA_SHAPE_PARAMETER,      self.alpha_parameter)
        self.spheres_model_part.ProcessInfo.SetValue(INTEGRAL_TOLERANCE,         self.integral_tolerance)

        # Interstitial fluid properties
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_DENSITY,              self.fluid_density)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_VISCOSITY,            self.fluid_viscosity)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_THERMAL_CONDUCTIVITY, self.fluid_thermal_conductivity)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_HEAT_CAPACITY,        self.fluid_heat_capacity)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_TEMPERATURE,          self.fluid_temperature)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_VELOCITY,             self.fluid_velocity)

    #----------------------------------------------------------------------------------------------
    def InitializeCPlusPlusUtilities(self):
        self.thermal_data_utils.ExecuteInitialize(self.spheres_model_part,self.fem_model_part)

        if (self.compute_voronoi or self.compute_porosity):
            self.tesselation_utils.ExecuteInitialize(self.spheres_model_part, self.compute_voronoi, self.compute_porosity)

        if (self.write_graph):
            self.graph_utils.ExecuteInitialize(self.PostGraphParticleTempMin,
                                               self.PostGraphParticleTempMax,
                                               self.PostGraphParticleTempAvg,
                                               self.PostGraphParticleTempDev,
                                               self.PostGraphModelTempAvg,
                                               self.PostGraphFluxContributions)

    #----------------------------------------------------------------------------------------------
    def IsTimeToUpdateVoronoi(self):
        if (self.compute_voronoi):
            step = self.spheres_model_part.ProcessInfo[TIME_STEPS]
            freq = self.voronoi_tesselation_frequency
            return step == 1 or (freq != 0 and step%freq == 0)
        else:
            return False

    #----------------------------------------------------------------------------------------------
    def IsTimeToUpdatePorosity(self):
        if (self.compute_porosity):
            step = self.spheres_model_part.ProcessInfo[TIME_STEPS]
            freq = self.porosity_update_frequency
            return step == 1 or (freq != 0 and step%freq == 0)
        else:
            return False
