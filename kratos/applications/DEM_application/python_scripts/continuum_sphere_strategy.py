from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import math

import sphere_strategy as SolverStrategy
BaseExplicitStrategy = SolverStrategy.ExplicitStrategy

class ExplicitStrategy(BaseExplicitStrategy):

    def __init__(self, model_part, fem_model_part, cluster_model_part, inlet_model_part, creator_destructor, dem_fem_search, scheme, Param, procedures):

        BaseExplicitStrategy.__init__(self, model_part, fem_model_part, cluster_model_part, inlet_model_part, creator_destructor, dem_fem_search, scheme, Param, procedures)

        self.print_skin_sphere = self.Var_Translator(Param.PostSkinSphere)

        if (self.delta_option > 0):
            self.case_option = 2     #MSIMSI. only 2 cases, with delta or without but continuum always.

        self.fixed_vel_top = Param.LoadingVelocityTop
        self.fixed_vel_bot = Param.LoadingVelocityBot

        if (self.Var_Translator(Param.DontSearchUntilFailure)):
            print ("Search is not active until a bond is broken.")
            self.search_control = 0
            if (len(fem_model_part.Nodes) > 0 or Param.TestType== "BTS"):   #MSI. This activates the search since there are fem contact elements. however only the particle - fem search should be active.
                print ("WARNING: Search should be activated since there might contact with FEM.")

        self.test_type = Param.TestType

        self.amplified_continuum_search_radius_extension = Param.AmplifiedSearchRadiusExtension

    def Initialize(self):
        self.SetVariablesAndOptions()

        # ADDITIONAL VARIABLES AND OPTIONS
        self.model_part.ProcessInfo.SetValue(AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION, self.amplified_continuum_search_radius_extension)
        self.model_part.ProcessInfo.SetValue(CONTACT_MESH_OPTION, self.contact_mesh_option)

        if ((self.test_type == "Triaxial") or (self.test_type == "Hydrostatic")):
            self.model_part.ProcessInfo.SetValue(TRIAXIAL_TEST_OPTION, 1)
        else:
            self.model_part.ProcessInfo.SetValue(TRIAXIAL_TEST_OPTION, 0)

        self.model_part.ProcessInfo.SetValue(FIXED_VEL_TOP, self.fixed_vel_top)
        self.model_part.ProcessInfo.SetValue(FIXED_VEL_BOT, self.fixed_vel_bot)
        ##################################

        if (self.Parameters.IntegrationScheme == 'Verlet_Velocity'):
            self.cplusplus_strategy = ContinuumVerletVelocitySolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                                            self.delta_option, self.creator_destructor, self.dem_fem_search, self.time_integration_scheme, self.search_strategy)
        else:
            self.cplusplus_strategy = ContinuumExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                  self.delta_option, self.creator_destructor, self.dem_fem_search, self.time_integration_scheme, self.search_strategy)

        self.cplusplus_strategy.Initialize()  # Calls the cplusplus_strategy Initialize function (initializes all elements and performs other necessary tasks before starting the time loop) (C++)

    def Initial_Critical_Time(self):
        (self.cplusplus_strategy).InitialTimeStepCalculation()

    def PrepareContactElementsForPrinting(self):
        (self.cplusplus_strategy).PrepareContactElementsForPrinting()

    def AddAdditionalVariables(self, model_part, Param):
        model_part.AddNodalSolutionStepVariable(COHESIVE_GROUP)  # Continuum group
        model_part.AddNodalSolutionStepVariable(SKIN_SPHERE)

    def ModifyProperties(self, properties):
        BaseExplicitStrategy.ModifyProperties(self, properties)

        ContinuumConstitutiveLawString = properties[DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME]
        ContinuumConstitutiveLaw = globals().get(ContinuumConstitutiveLawString)()
        ContinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties)
