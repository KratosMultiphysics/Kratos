from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import sphere_strategy as SolverStrategy
BaseExplicitStrategy = SolverStrategy.ExplicitStrategy

import math

class ExplicitStrategy(BaseExplicitStrategy):

    def __init__(self, spheres_model_part, fem_model_part, cluster_model_part, inlet_model_part, creator_destructor, dem_fem_search, scheme, Param, procedures):

        BaseExplicitStrategy.__init__(self, spheres_model_part, fem_model_part, cluster_model_part, inlet_model_part, creator_destructor, dem_fem_search, scheme, Param, procedures)

        self.print_skin_sphere = self.Var_Translator(Param.PostSkinSphere)

        if (self.delta_option > 0):
            self.case_option = 2     #MSIMSI. only 2 cases, with delta or without but continuum always.

        if not hasattr(Param, "LoadingVelocityTop"):
            self.fixed_vel_top = 0
        else:
            self.fixed_vel_top = Param.LoadingVelocityTop

        if not hasattr(Param, "LoadingVelocityBot"):
            self.fixed_vel_bot = 0
        else:
            self.fixed_vel_bot = Param.LoadingVelocityBot

        if (self.Var_Translator(Param.DontSearchUntilFailure)):
            print ("Search is not active until a bond is broken.")
            self.search_control = 0
            if (len(fem_model_part.Nodes) > 0 or Param.TestType== "BTS"):   #MSI. This activates the search since there are fem contact elements. however only the particle - fem search should be active.
                print ("WARNING: Search should be activated since there might contact with FEM.")

        if not hasattr(Param, "TestType"):
            self.test_type = "None"
        else:
            self.test_type = Param.TestType

        self.amplified_continuum_search_radius_extension = Param.AmplifiedSearchRadiusExtension
        
        if hasattr(Param, "MaxAmplificationRatioOfSearchRadius"):
            self.max_amplification_ratio_of_search_radius = Param.MaxAmplificationRatioOfSearchRadius
        else:
            self.max_amplification_ratio_of_search_radius = 0.0
            
        if not hasattr(Param, "PostPoissonRatio"):
            self.poisson_ratio_option = 0
        else:
            self.poisson_ratio_option = self.Var_Translator(Param.PostPoissonRatio)
            
        if not hasattr(Param, "PoissonEffectOption"):
            self.poisson_effect_option = 0
        else:
            self.poisson_effect_option = self.Var_Translator(Param.PoissonEffectOption)

        if not hasattr(Param, "ShearStrainParallelToBondOption"):
            self.shear_strain_parallel_to_bond_option = 0
        else:
            self.shear_strain_parallel_to_bond_option = self.Var_Translator(Param.ShearStrainParallelToBondOption)

        if (self.poisson_effect_option or self.shear_strain_parallel_to_bond_option):
            self.compute_stress_tensor_option = 1

    def CreateCPlusPlusStrategy(self):
        self.SetVariablesAndOptions()

        # ADDITIONAL VARIABLES AND OPTIONS
        self.spheres_model_part.ProcessInfo.SetValue(AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION, self.amplified_continuum_search_radius_extension)
        self.spheres_model_part.ProcessInfo.SetValue(MAX_AMPLIFICATION_RATIO_OF_THE_SEARCH_RADIUS, self.max_amplification_ratio_of_search_radius)
        self.spheres_model_part.ProcessInfo.SetValue(CONTACT_MESH_OPTION, self.contact_mesh_option)
        self.spheres_model_part.ProcessInfo.SetValue(COMPUTE_STRESS_TENSOR_OPTION, self.compute_stress_tensor_option)

        if ((self.test_type == "Triaxial") or (self.test_type == "Hydrostatic")):
            self.spheres_model_part.ProcessInfo.SetValue(TRIAXIAL_TEST_OPTION, 1)
        else:
            self.spheres_model_part.ProcessInfo.SetValue(TRIAXIAL_TEST_OPTION, 0)

        self.spheres_model_part.ProcessInfo.SetValue(FIXED_VEL_TOP, self.fixed_vel_top)
        self.spheres_model_part.ProcessInfo.SetValue(FIXED_VEL_BOT, self.fixed_vel_bot)
        
        self.spheres_model_part.ProcessInfo.SetValue(POISSON_EFFECT_OPTION, self.poisson_effect_option)
        self.spheres_model_part.ProcessInfo.SetValue(SHEAR_STRAIN_PARALLEL_TO_BOND_OPTION, self.shear_strain_parallel_to_bond_option)

        ##################################

        if (self.Parameters.IntegrationScheme == 'Verlet_Velocity'):
            self.cplusplus_strategy = ContinuumVerletVelocitySolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                                            self.delta_option, self.creator_destructor, self.dem_fem_search, self.time_integration_scheme, self.search_strategy)
        else:
            self.cplusplus_strategy = ContinuumExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                  self.delta_option, self.creator_destructor, self.dem_fem_search, self.time_integration_scheme, self.search_strategy)
    def Initialize(self):        
        self.cplusplus_strategy.Initialize()  # Calls the cplusplus_strategy Initialize function (initializes all elements and performs other necessary tasks before starting the time loop) (C++)

    def SetContinuumType(self):
        self.continuum_type = True
    
    def Initial_Critical_Time(self):
        (self.cplusplus_strategy).InitialTimeStepCalculation()

    def PrepareContactElementsForPrinting(self):
        (self.cplusplus_strategy).PrepareContactElementsForPrinting()

    def AddAdditionalVariables(self, spheres_model_part, Param):
        spheres_model_part.AddNodalSolutionStepVariable(COHESIVE_GROUP)  # Continuum group
        spheres_model_part.AddNodalSolutionStepVariable(SKIN_SPHERE)

    def ModifyProperties(self, properties):
        BaseExplicitStrategy.ModifyProperties(self, properties)

        ContinuumConstitutiveLawString = properties[DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME]
        ContinuumConstitutiveLaw = globals().get(ContinuumConstitutiveLawString)()
        ContinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties)
