from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *


# DEM Application using Velocity Verlet 2-step scheme for continuum

import continuum_sphere_strategy as SolverStrategy

BaseExplicitStrategy = SolverStrategy.ExplicitStrategy

class ExplicitStrategy(BaseExplicitStrategy):

    def __init__(self, model_part, fem_model_part, cluster_model_part, inlet_model_part, creator_destructor, dem_fem_search, DEM_parameters, procedures):

        BaseExplicitStrategy.__init__(self, model_part, fem_model_part, cluster_model_part, inlet_model_part, creator_destructor, dem_fem_search, DEM_parameters, procedures)

    def AddAdditionalVariables(self, model_part, DEM_parameters):

        BaseExplicitStrategy.AddAdditionalVariables(self, model_part, DEM_parameters)

    def CreateCPlusPlusStrategy(self):

        #BaseExplicitStrategy.Initialize  (revisar si es pot cridar desde el basetype)
        #self.cplusplus_strategy = IterativeExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,

        # Setting ProcessInfo variables

        # SIMULATION FLAGS
        self.model_part.ProcessInfo.SetValue(VIRTUAL_MASS_OPTION, self.virtual_mass_option)
        self.model_part.ProcessInfo.SetValue(CRITICAL_TIME_OPTION, self.critical_time_option)
        self.model_part.ProcessInfo.SetValue(CASE_OPTION, self.case_option)
        self.model_part.ProcessInfo.SetValue(TRIHEDRON_OPTION, self.trihedron_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, ROTATION_OPTION, self.rotation_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, BOUNDING_BOX_OPTION, self.bounding_box_option) #TODO: check that this is finding the function in base classes
        self.model_part.ProcessInfo.SetValue(SEARCH_CONTROL, self.search_control)
        self.model_part.ProcessInfo.SetValue(FIX_VELOCITIES_FLAG, self.fix_velocities_flag)
        self.model_part.ProcessInfo.SetValue(NEIGH_INITIALIZED, 0)
        self.model_part.ProcessInfo.SetValue(CLEAN_INDENT_OPTION, self.clean_init_indentation_option)
        self.model_part.ProcessInfo.SetValue(COMPUTE_STRESS_TENSOR_OPTION, self.compute_stress_tensor_option)

        # GLOBAL PHISICAL ASPECTS
        self.model_part.ProcessInfo.SetValue(GRAVITY, self.gravity)

        # GLOBAL MATERIAL PROPERTIES
        self.model_part.ProcessInfo.SetValue(NODAL_MASS_COEFF, self.nodal_mass_coeff)

        # PRINTING VARIABLES
        self.model_part.ProcessInfo.SetValue(PRINT_EXPORT_ID, self.print_export_id)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, ROLLING_FRICTION_OPTION, self.rolling_friction_option)

        # TIME RELATED PARAMETERS
        self.model_part.ProcessInfo.SetValue(DELTA_TIME, self.delta_time)

        for properties in self.model_part.Properties:
            self.ModifyProperties(properties)

        for properties in self.inlet_model_part.Properties:
            self.ModifyProperties(properties)

        # CONTINUUM
        self.search_increment_for_walls = 0.0 # for the moment, until all bugs are removed
        self.model_part.ProcessInfo.SetValue(SEARCH_RADIUS_INCREMENT, self.search_increment)
        self.spheres_model_part.ProcessInfo.SetValue(SEARCH_RADIUS_INCREMENT_FOR_WALLS, self.search_increment_for_walls)
        self.model_part.ProcessInfo.SetValue(AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION, self.amplified_continuum_search_radius_extension)
        self.model_part.ProcessInfo.SetValue(LOCAL_RESOLUTION_METHOD, self.local_resolution_method)

        if self.contact_mesh_option:
            self.model_part.ProcessInfo.SetValue(CONTACT_MESH_OPTION, 1) #TODO: convert this variable to BOOL in Kratos
        else:
            self.model_part.ProcessInfo.SetValue(CONTACT_MESH_OPTION, 0)

        #self.model_part.ProcessInfo.SetValue(FAILURE_CRITERION_OPTION, self.failure_criterion_option)

        if ( (self.test_type == "Triaxial") or (self.test_type == "Hydrostatic")):
            self.model_part.ProcessInfo.SetValue(TRIAXIAL_TEST_OPTION, 1)

        # RESOLUTION METHODS AND PARAMETERS
        # Creating the solution strategy
        self.settings = ExplicitSolverSettings()
        self.settings.r_model_part = self.model_part
        self.settings.contact_model_part = self.contact_model_part
        self.settings.fem_model_part = self.fem_model_part
        self.settings.inlet_model_part = self.inlet_model_part
        self.settings.cluster_model_part = self.cluster_model_part

        self.cplusplus_strategy = VelocityVerletSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor, self.delta_option, self.creator_destructor, self.dem_fem_search, self.search_strategy)

    def Initialize(self):

        self.cplusplus_strategy.Initialize()  # Calls the cplusplus_strategy Initialize function (initializes all elements and performs other necessary tasks before iterating)
