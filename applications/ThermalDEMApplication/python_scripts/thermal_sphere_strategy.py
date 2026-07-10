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
    def AddVariables(self):
        # Add standard variables
        BaseStrategy.AddVariables(self)

        # Add thermal variables to all model parts
        self.AddThermalVariables()
        
    #----------------------------------------------------------------------------------------------
    def CreateCPlusPlusStrategy(self):
        # Set standard options
        BaseStrategy.SetVariablesAndOptions(self)

        # Set thermal options (set ProcessInfo values)
        self.SetThermalVariablesAndOptions()

        # Create cpp strategy object
        self.CreateCPlusPlusThermalStrategy()

    #----------------------------------------------------------------------------------------------
    def ModifyProperties(self, properties, param = 0):
        if param:
            return

        # Set standard properties
        BaseStrategy.ModifyProperties(self, properties, param)

        # Set pointers: constitutive laws (heat transfer models) / time integration scheme / numerical integration method / heat map utilities
        self.SetConstitutiveLaw(properties)
        self.SetThermalIntegrationScheme(properties)
        self.SetNumericalIntegrationMethod(properties)
    
    #----------------------------------------------------------------------------------------------
    def Initialize(self):
        # Initialize utilities
        # (important to be before the initialization of elements, because temperature is set here)
        self.InitializeCPlusPlusUtilities()

        # Base class initializer
        # (initialize the strategy and the elements, so temperature must be already set at this point)
        BaseStrategy.Initialize(self)
    
    #----------------------------------------------------------------------------------------------
    def InitializeSolutionStep(self):
        if (self.compute_forces_option):
            BaseStrategy.InitializeSolutionStep(self)
        else:
            (self.cplusplus_strategy).InitializeSolutionStep()

        # Perform tesselation-dependent tasks (triangulation or tetrahedralization)
        if (self.IsTimeToUpdateVoronoi() or self.IsTimeToUpdatePorosity()):
            self.tesselation_utils.ExecuteInitializeSolutionStep(self.spheres_model_part)
    
    #----------------------------------------------------------------------------------------------
    def Predict(self):
        if (self.compute_forces_option):
            BaseStrategy.Predict(self)
    
    #----------------------------------------------------------------------------------------------
    def SolveSolutionStep(self):
        # Solve step according to motion type
        if (self.compute_forces_option):
            (self.cplusplus_strategy).SolveSolutionStep()
        else:
            (self.cplusplus_strategy).SolveSolutionStepStatic()
        
        return True
    
    #----------------------------------------------------------------------------------------------
    def FinalizeSolutionStep(self):
        BaseStrategy.FinalizeSolutionStep(self)
        
        # Write output graphs
        if (self.write_graph):
            self.graph_utils.ExecuteFinalizeSolutionStep(self.spheres_model_part)

        # Merge particle heat maps to global heat maps
        if (self.PostHeatMapGeneration):
            self.heat_map_utils.ExecuteFinalizeSolutionStep(self.spheres_model_part)

    #----------------------------------------------------------------------------------------------
    def Finalize(self):
        BaseStrategy.Finalize(self)

        # Close graph files
        if (self.write_graph):
            self.graph_utils.ExecuteFinalize()

        # Write global heat maps
        if (self.PostHeatMapGeneration):
            self.heat_map_utils.ExecuteFinalize(self.spheres_model_part)

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
        self.compute_forces_option       = self.thermal_settings["compute_forces"].GetBool()
        self.compute_motion_option       = self.thermal_settings["compute_motion"].GetBool()
        self.auto_solve_frequency_option = self.thermal_settings["automatic_solve_frequency"].GetBool()

        # Frequencies
        self.thermal_solve_frequency       = self.thermal_settings["thermal_solve_frequency"].GetInt()
        self.voronoi_tesselation_frequency = self.thermal_settings["voronoi_tesselation_frequency"].GetInt()
        self.porosity_update_frequency     = self.thermal_settings["porosity_update_frequency"].GetInt()

        # Integration scheme and method
        self.thermal_integration_scheme   = self.thermal_settings["thermal_integration_scheme"].GetString()
        self.numerical_integration_method = self.thermal_settings["numerical_integration_method"].GetString()

        # Models for heat transfer
        self.direct_conduction_model   = self.thermal_settings["direct_conduction_model"].GetString()
        self.indirect_conduction_model = self.thermal_settings["indirect_conduction_model"].GetString()
        self.nusselt_correlation       = self.thermal_settings["nusselt_correlation"].GetString()
        self.radiation_model           = self.thermal_settings["radiation_model"].GetString()
        self.adjusted_contact_model    = self.thermal_settings["adjusted_contact_model"].GetString()
        self.voronoi_method            = self.thermal_settings["voronoi_method"].GetString()
        self.porosity_method           = self.thermal_settings["porosity_method"].GetString()

        self.heat_generation_model = []
        for model in self.thermal_settings["heat_generation_model"].values():
            self.heat_generation_model.append(model.GetString())

        # Active heat transfer mechanisms
        self.compute_direct_conduction_option   = GetBoolParameterIfItExists(self.thermal_settings, "compute_direct_conduction")
        self.compute_indirect_conduction_option = GetBoolParameterIfItExists(self.thermal_settings, "compute_indirect_conduction")
        self.compute_convection_option          = GetBoolParameterIfItExists(self.thermal_settings, "compute_convection")
        self.compute_radiation_option           = GetBoolParameterIfItExists(self.thermal_settings, "compute_radiation")
        self.compute_heat_generation_option     = GetBoolParameterIfItExists(self.thermal_settings, "compute_heat_generation")
        self.compute_adjusted_contact_option    = GetBoolParameterIfItExists(self.thermal_settings, "compute_adjusted_contact")
        
        # Model parameters
        self.min_conduction_distance  = self.thermal_settings["min_conduction_distance"].GetDouble()
        self.max_conduction_distance  = self.thermal_settings["max_conduction_distance"].GetDouble()
        self.conduction_radius        = self.thermal_settings["conduction_radius"].GetDouble()
        self.fluid_layer_thickness    = self.thermal_settings["fluid_layer_thickness"].GetDouble()
        self.isothermal_core_radius   = self.thermal_settings["isothermal_core_radius"].GetDouble()
        self.max_radiation_distance   = self.thermal_settings["max_radiation_distance"].GetDouble()
        self.heat_generation_ratio    = self.thermal_settings["heat_generation_ratio"].GetDouble()
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
        
        # Post options
        self.PostGraphParticleTempAll             = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphParticleTempAll")
        self.PostGraphParticleTempMin             = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphParticleTempMin")
        self.PostGraphParticleTempMax             = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphParticleTempMax")
        self.PostGraphParticleTempAvg             = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphParticleTempAvg")
        self.PostGraphParticleTempAvgVol          = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphParticleTempAvgVol")
        self.PostGraphParticleTempDev             = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphParticleTempDev")
        self.PostGraphMechanicalEnergy            = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphMechanicalEnergy")
        self.PostGraphDissipatedEnergy            = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphDissipatedEnergy")
        self.PostGraphThermalEnergy               = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphThermalEnergy")
        self.PostGraphHeatFluxContributions       = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphHeatFluxContributions")
        self.PostGraphHeatGenerationValues        = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphHeatGenerationValues")
        self.PostGraphHeatGenerationContributions = GetBoolParameterIfItExists(self.DEM_parameters, "PostGraphHeatGenerationContributions")
        self.PostHeatMapGeneration                = GetBoolParameterIfItExists(self.DEM_parameters, "PostHeatMapGeneration")

        self.heat_map_corner1         = Vector(3)
        self.heat_map_corner1[0]      = min(self.thermal_settings["heat_map_corners"][0][0].GetDouble(),self.thermal_settings["heat_map_corners"][1][0].GetDouble())
        self.heat_map_corner1[1]      = min(self.thermal_settings["heat_map_corners"][0][1].GetDouble(),self.thermal_settings["heat_map_corners"][1][1].GetDouble())
        self.heat_map_corner1[2]      = min(self.thermal_settings["heat_map_corners"][0][2].GetDouble(),self.thermal_settings["heat_map_corners"][1][2].GetDouble())
        self.heat_map_corner2         = Vector(3)
        self.heat_map_corner2[0]      = max(self.thermal_settings["heat_map_corners"][0][0].GetDouble(),self.thermal_settings["heat_map_corners"][1][0].GetDouble())
        self.heat_map_corner2[1]      = max(self.thermal_settings["heat_map_corners"][0][1].GetDouble(),self.thermal_settings["heat_map_corners"][1][1].GetDouble())
        self.heat_map_corner2[2]      = max(self.thermal_settings["heat_map_corners"][0][2].GetDouble(),self.thermal_settings["heat_map_corners"][1][2].GetDouble())
        self.heat_map_subdivisions    = Vector(3)
        self.heat_map_subdivisions[0] = self.thermal_settings["heat_map_subdivisions"][0].GetInt()
        self.heat_map_subdivisions[1] = self.thermal_settings["heat_map_subdivisions"][1].GetInt()
        self.heat_map_subdivisions[2] = self.thermal_settings["heat_map_subdivisions"][2].GetInt()

    #----------------------------------------------------------------------------------------------
    def CheckProjectParameters(self):
        # Forces and motion calculation
        if (self.compute_motion_option == True and self.compute_forces_option == False):
            Logger.PrintWarning('ThermalDEM', '\nActivation of "compute_motion" requires "compute_forces" to be enabled. The simulation will run with "compute_forces" enabled.\n')
        
        # Time integration scheme
        if (self.thermal_integration_scheme != "forward_euler"):
            raise Exception('ThermalDEM', 'Time integration scheme \'' + self.thermal_integration_scheme + '\' is not implemented.') 
        
        # Numerical integration method
        if (self.numerical_integration_method != "adaptive_simpson"):
            raise Exception('ThermalDEM', 'Numerical integration method \'' + self.numerical_integration_method + '\' is not implemented.')
        
        # Heat transfer models
        if (self.direct_conduction_model != "batchelor_obrien_simple"   and
            self.direct_conduction_model != "batchelor_obrien_complete" and
            self.direct_conduction_model != "batchelor_obrien_modified" and
            self.direct_conduction_model != "thermal_pipe"              and
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

        for model in self.heat_generation_model:
            if (model != "sliding_friction" and
                model != "rolling_friction" and
                model != "contact_damping"):
                raise Exception('ThermalDEM', 'Heat generation model \'' + model + '\' is not implemented.')

        if (self.adjusted_contact_model != "zhou"             and
            self.adjusted_contact_model != "lu"               and
            self.adjusted_contact_model != "morris_area"      and
            self.adjusted_contact_model != "morris_area_time" and
            self.adjusted_contact_model != "rangel_area"      and
            self.adjusted_contact_model != "rangel_area_time"):
            raise Exception('ThermalDEM', 'Adjusted contact model \'' + self.adjusted_contact_model + '\' is not implemented.')

        # Other methods
        if (self.voronoi_method != "tesselation" and
            self.voronoi_method != "porosity"):
            raise Exception('ThermalDEM', 'Voronoi method \'' + self.voronoi_method + '\' is not implemented.')
        
        if (self.porosity_method != "global"              and
            self.porosity_method != "average_convex_hull" and
            self.porosity_method != "average_alpha_shape"):
            raise Exception('ThermalDEM', 'Porosity method \'' + self.porosity_method + '\' is not implemented.')

        # Model parameters values
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
        if (self.conduction_radius < 0):
            self.conduction_radius = 0
        if (self.fluid_layer_thickness < 0):
            self.fluid_layer_thickness = 0
        if (self.isothermal_core_radius < 0):
            self.isothermal_core_radius = 0
        if (self.isothermal_core_radius > 1):
            self.isothermal_core_radius = 1
        if (self.max_radiation_distance < 0 ):
            self.max_radiation_distance = 0
        if (self.global_porosity < 0 or self.global_porosity >= 1):
            raise Exception('ThermalDEM', '"global_porosity" must be between zero and one.')
        if (self.alpha_parameter < 0):
            raise Exception('ThermalDEM', '"alpha_shape_parameter" must be positive.')
        if (self.integral_tolerance <= 0):
            raise Exception('ThermalDEM', '"integral_tolerance" must be positive.')
        
        # Fluid properties values
        if (self.fluid_density              <= 0 or
            self.fluid_viscosity            <= 0 or
            self.fluid_thermal_conductivity <= 0 or
            self.fluid_heat_capacity        <= 0):
            raise Exception('ThermalDEM', '"global_fluid_properties" must contain positive values for material properties.')
        
        # Post options
        if (self.heat_map_corner1[0] == self.heat_map_corner2[0] or
            self.heat_map_corner1[1] == self.heat_map_corner2[1] or
            self.heat_map_corner1[2] == self.heat_map_corner2[2]):
            raise Exception('ThermalDEM', '"heat_map_corners" must contain two vectors with the X,Y,Z coordinates of points that define opposite corners of a cuboid.')
        
        if (self.heat_map_subdivisions[0] < 1 or
            self.heat_map_subdivisions[1] < 1 or
            self.heat_map_subdivisions[2] < 1):
            raise Exception('ThermalDEM', '"heat_map_subdivisions" must contain a vector with 3 positive values for the number of subdivisions in X,Y,Z directions.')

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
              self.voronoi_method            == "porosity"    and
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
        if (self.PostGraphParticleTempAll       or
            self.PostGraphParticleTempMin       or
            self.PostGraphParticleTempMax       or
            self.PostGraphParticleTempAvg       or 
            self.PostGraphParticleTempAvgVol    or
            self.PostGraphParticleTempDev       or
            self.PostGraphMechanicalEnergy      or
            self.PostGraphDissipatedEnergy      or
            self.PostGraphThermalEnergy         or
            self.PostGraphHeatFluxContributions or
            self.PostGraphHeatGenerationValues  or
            self.PostGraphHeatGenerationContributions):
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

        if (self.PostHeatMapGeneration):
            self.heat_map_utils = HeatMapUtilities()
    
    #----------------------------------------------------------------------------------------------
    def AddThermalVariables(self):
        self.spheres_model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        self.cluster_model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        self.inlet_model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        self.fem_model_part.AddNodalSolutionStepVariable(TEMPERATURE)

        self.spheres_model_part.AddNodalSolutionStepVariable(HEATFLUX)
        self.cluster_model_part.AddNodalSolutionStepVariable(HEATFLUX)
        self.inlet_model_part.AddNodalSolutionStepVariable(HEATFLUX)
        self.fem_model_part.AddNodalSolutionStepVariable(HEATFLUX)

    #----------------------------------------------------------------------------------------------
    def SetConstitutiveLaw(self, properties):
        # Direct conduction
        if self.direct_conduction_model == "batchelor_obrien_simple":
            class_name = "DirectConductionBOBSimple"
        elif self.direct_conduction_model == "batchelor_obrien_complete":
            class_name = "DirectConductionBOBComplete"
        elif self.direct_conduction_model == "batchelor_obrien_modified":
            class_name = "DirectConductionBOBModified"
        elif self.direct_conduction_model == "thermal_pipe":
            class_name = "DirectConductionPipe"
        elif self.direct_conduction_model == "collisional":
            class_name = "DirectConductionCollision"
        else:
            raise Exception('ThermalDEM', 'Direct thermal conduction model \'' + self.direct_conduction_model + '\' is not implemented.')
        try:
            object = eval(class_name)()
        except:
            raise Exception('The class corresponding to the direct thermal conduction model named ' + class_name + ' has not been added to python. Please, select a different name or add the required class.')
        object.SetHeatExchangeMechanismInProperties(properties, True)

        # Indirect conduction
        if self.indirect_conduction_model == "surrounding_layer":
            class_name = "IndirectConductionSurroundLayer"
        elif self.indirect_conduction_model == "voronoi_a":
            class_name = "IndirectConductionVoronoiA"
        elif self.indirect_conduction_model == "voronoi_b":
            class_name = "IndirectConductionVoronoiB"
        elif self.indirect_conduction_model == "vargas_mccarthy":
            class_name = "IndirectConductionVargas"
        else:
            raise Exception('ThermalDEM', 'Indirect thermal conduction model \'' + self.indirect_conduction_model + '\' is not implemented.')
        try:
            object = eval(class_name)()
        except:
            raise Exception('The class corresponding to the indirect thermal conduction model named ' + class_name + ' has not been added to python. Please, select a different name or add the required class.')
        object.SetHeatExchangeMechanismInProperties(properties, True)

        # Convection
        if self.nusselt_correlation == "sphere_hanz_marshall":
            class_name = "NusseltHanzMarshall"
        elif self.nusselt_correlation == "sphere_whitaker":
            class_name = "NusseltWhitaker"
        elif self.nusselt_correlation == "sphere_gunn":
            class_name = "NusseltGunn"
        elif self.nusselt_correlation == "sphere_li_mason":
            class_name = "NusseltLiMason"
        else:
            raise Exception('ThermalDEM', 'Nusselt number correlation \'' + self.nusselt_correlation + '\' is not implemented.')
        try:
            object = eval(class_name)()
        except:
            raise Exception('The class corresponding to the nusselt number correlation named ' + class_name + ' has not been added to python. Please, select a different name or add the required class.')
        object.SetHeatExchangeMechanismInProperties(properties, True)

        # Radiation
        if self.radiation_model == "continuum_zhou":
            class_name = "RadiationContinuumZhou"
        elif self.radiation_model == "continuum_krause":
            class_name = "RadiationContinuumKrause"
        else:
            raise Exception('ThermalDEM', 'Thermal radiation model \'' + self.radiation_model + '\' is not implemented.')
        try:
            object = eval(class_name)()
        except:
            raise Exception('The class corresponding to the thermal radiation model named ' + class_name + ' has not been added to python. Please, select a different name or add the required class.')
        object.SetHeatExchangeMechanismInProperties(properties, True)

        # Heat generation
        class_name = "GenerationDissipation"
        try:
            object = eval(class_name)()
        except:
            raise Exception('The class corresponding to the heat generation model named ' + class_name + ' has not been added to python. Please, select a different name or add the required class.')
        object.SetHeatGenerationMechanismInProperties(properties, True)

        # Real contact
        if self.adjusted_contact_model == "zhou":
            class_name = "RealContactZhou"
        elif self.adjusted_contact_model == "lu":
            class_name = "RealContactLu"
        elif self.adjusted_contact_model == "morris_area":
            class_name = "RealContactMorrisArea"
        elif self.adjusted_contact_model == "morris_area_time":
            class_name = "RealContactMorrisAreaTime"
        elif self.adjusted_contact_model == "rangel_area":
            class_name = "RealContactRangelArea"
        elif self.adjusted_contact_model == "rangel_area_time":
            class_name = "RealContactRangelAreaTime"
        else:
            raise Exception('ThermalDEM', 'Real contact model \'' + self.adjusted_contact_model + '\' is not implemented.')
        try:
            object = eval(class_name)()
        except:
            raise Exception('The class corresponding to the real contact model named ' + class_name + ' has not been added to python. Please, select a different name or add the required class.')
        object.SetRealContactModelInProperties(properties, True)

    #----------------------------------------------------------------------------------------------
    def SetThermalIntegrationScheme(self, properties):
        if properties.Has(THERMAL_INTEGRATION_SCHEME_NAME):
            input_name = properties[THERMAL_INTEGRATION_SCHEME_NAME]
        else:
            input_name = self.thermal_integration_scheme

        if input_name == "forward_euler":
            class_name = "ThermalForwardEulerScheme"
        else:
            raise Exception('ThermalDEM', 'Time integration scheme \'' + input_name + '\' is not implemented.')

        try:
            object = eval(class_name)()
        except:
            raise Exception('The class corresponding to the time integration scheme named ' + class_name + ' has not been added to python. Please, select a different name or add the required class.')
        
        object.SetThermalIntegrationSchemeInProperties(properties, True)

    #----------------------------------------------------------------------------------------------
    def SetNumericalIntegrationMethod(self, properties):
        if properties.Has(NUMERICAL_INTEGRATION_METHOD_NAME):
            input_name = properties[NUMERICAL_INTEGRATION_METHOD_NAME]
        else:
            input_name = self.numerical_integration_method

        if input_name == "adaptive_simpson":
            class_name = "AdaptiveSimpsonQuadrature"
        else:
            raise Exception('ThermalDEM', 'Numerical integration method \'' + input_name + '\' is not implemented.')

        try:
            object = eval(class_name)()
        except:
            raise Exception('The class corresponding to the numerical integration method named ' + class_name + ' has not been added to python. Please, select a different name or add the required class.')
        
        object.SetNumericalIntegrationMethodInProperties(properties, True)

    #----------------------------------------------------------------------------------------------
    def SetThermalVariablesAndOptions(self):
        # General options
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, COMPUTE_FORCES_OPTION, self.compute_forces_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, COMPUTE_MOTION_OPTION, self.compute_motion_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, AUTO_SOLVE_FREQUENCY_OPTION, self.auto_solve_frequency_option)
        self.spheres_model_part.ProcessInfo.SetValue(THERMAL_FREQUENCY, self.thermal_solve_frequency)

        # Models for heat transfer
        self.spheres_model_part.ProcessInfo.SetValue(DIRECT_CONDUCTION_MODEL_NAME,   self.direct_conduction_model)
        self.spheres_model_part.ProcessInfo.SetValue(INDIRECT_CONDUCTION_MODEL_NAME, self.indirect_conduction_model)
        self.spheres_model_part.ProcessInfo.SetValue(CONVECTION_MODEL_NAME,          self.nusselt_correlation)
        self.spheres_model_part.ProcessInfo.SetValue(RADIATION_MODEL_NAME,           self.radiation_model)
        self.spheres_model_part.ProcessInfo.SetValue(REAL_CONTACT_MODEL_NAME,        self.adjusted_contact_model)
        self.spheres_model_part.ProcessInfo.SetValue(VORONOI_METHOD_NAME,            self.voronoi_method)
        self.spheres_model_part.ProcessInfo.SetValue(POROSITY_METHOD_NAME,           self.porosity_method)

        # Active heat transfer mechanisms
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, DIRECT_CONDUCTION_OPTION,   self.compute_direct_conduction_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, INDIRECT_CONDUCTION_OPTION, self.compute_indirect_conduction_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, CONVECTION_OPTION,          self.compute_convection_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, RADIATION_OPTION,           self.compute_radiation_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, HEAT_GENERATION_OPTION,     self.compute_heat_generation_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, GENERATION_SLIDING_OPTION,  self.compute_heat_generation_option and "sliding_friction" in self.heat_generation_model)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, GENERATION_ROLLING_OPTION,  self.compute_heat_generation_option and "rolling_friction" in self.heat_generation_model)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, GENERATION_DAMPING_OPTION,  self.compute_heat_generation_option and "contact_damping"  in self.heat_generation_model)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, REAL_CONTACT_OPTION,        self.compute_adjusted_contact_option)

        # Model parameters
        self.spheres_model_part.ProcessInfo.SetValue(MIN_CONDUCTION_DISTANCE,  self.min_conduction_distance)
        self.spheres_model_part.ProcessInfo.SetValue(MAX_CONDUCTION_DISTANCE,  self.max_conduction_distance)
        self.spheres_model_part.ProcessInfo.SetValue(CONDUCTION_RADIUS,        self.conduction_radius)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_LAYER_THICKNESS,    self.fluid_layer_thickness)
        self.spheres_model_part.ProcessInfo.SetValue(ISOTHERMAL_CORE_RADIUS,   self.isothermal_core_radius)
        self.spheres_model_part.ProcessInfo.SetValue(MAX_RADIATION_DISTANCE,   self.max_radiation_distance)
        self.spheres_model_part.ProcessInfo.SetValue(HEAT_GENERATION_RATIO,    self.heat_generation_ratio)
        self.spheres_model_part.ProcessInfo.SetValue(AVERAGE_POROSITY,         self.global_porosity)
        self.spheres_model_part.ProcessInfo.SetValue(ALPHA_SHAPE_PARAMETER,    self.alpha_parameter)
        self.spheres_model_part.ProcessInfo.SetValue(INTEGRAL_TOLERANCE,       self.integral_tolerance)

        # Interstitial fluid properties
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_DENSITY,              self.fluid_density)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_VISCOSITY,            self.fluid_viscosity)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_THERMAL_CONDUCTIVITY, self.fluid_thermal_conductivity)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_HEAT_CAPACITY,        self.fluid_heat_capacity)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_TEMPERATURE,          self.fluid_temperature)
        self.spheres_model_part.ProcessInfo.SetValue(FLUID_VELOCITY,             self.fluid_velocity)

        # Post options
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, HEAT_MAP_GENERATION_OPTION, self.PostHeatMapGeneration)
        self.spheres_model_part.ProcessInfo.SetValue(HEAT_MAP_COORDINATES_1, self.heat_map_corner1)
        self.spheres_model_part.ProcessInfo.SetValue(HEAT_MAP_COORDINATES_2, self.heat_map_corner2)
        self.spheres_model_part.ProcessInfo.SetValue(HEAT_MAP_SUBDIVISIONS,  self.heat_map_subdivisions)

    #----------------------------------------------------------------------------------------------
    def CreateCPlusPlusThermalStrategy(self):
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
    def InitializeCPlusPlusUtilities(self):
        self.thermal_data_utils.ExecuteInitialize(self.spheres_model_part,self.fem_model_part)

        if (self.compute_voronoi or self.compute_porosity):
            self.tesselation_utils.ExecuteInitialize(self.spheres_model_part, self.compute_voronoi, self.compute_porosity)

        if (self.write_graph):
            self.graph_utils.ExecuteInitialize(self.PostGraphParticleTempAll,
                                               self.PostGraphParticleTempMin,
                                               self.PostGraphParticleTempMax,
                                               self.PostGraphParticleTempAvg,
                                               self.PostGraphParticleTempAvgVol,
                                               self.PostGraphParticleTempDev,
                                               self.PostGraphMechanicalEnergy,
                                               self.PostGraphDissipatedEnergy,
                                               self.PostGraphThermalEnergy,
                                               self.PostGraphHeatFluxContributions,
                                               self.PostGraphHeatGenerationValues,
                                               self.PostGraphHeatGenerationContributions)

        if (self.PostHeatMapGeneration):
            self.heat_map_utils.ExecuteInitialize(self.spheres_model_part)

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
