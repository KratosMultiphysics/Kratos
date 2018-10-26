from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import sphere_strategy as SolverStrategy
BaseExplicitStrategy = SolverStrategy.ExplicitStrategy

import math

class ExplicitStrategy(BaseExplicitStrategy):

    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures):

        BaseExplicitStrategy.__init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures)

        self.print_skin_sphere = 0 #TODO: check if this variable is important. There's a similar one in DEM_procedures called PostSkinSphere
        if "PostSkinSphere" in DEM_parameters.keys():
            self.print_skin_sphere = DEM_parameters["PostSkinSphere"].GetBool()

        if (self.delta_option > 0):
            self.case_option = 2     #MSIMSI. only 2 cases, with delta or without but continuum always.

        if "DontSearchUntilFailure" in DEM_parameters.keys(): #TODO: important Todo. When Json gets divided in encapsulated parts, all these checks should be done in one functions, comparing with defaults!
            if DEM_parameters["DontSearchUntilFailure"].GetBool():
                print ("Search is not active until a bond is broken.")
                self.search_control = 0
                if (len(fem_model_part.Nodes) > 0 or DEM_parameters["TestType"].GetString() == "BTS"):   #MSI. This activates the search since there are fem contact elements. however only the particle - fem search should be active.
                    print ("WARNING: Search should be activated since there might contact with FEM.")

        if not "TestType" in DEM_parameters.keys():
            self.test_type = "None"
        else:
            self.test_type = DEM_parameters["TestType"].GetString()

        self.amplified_continuum_search_radius_extension = DEM_parameters["AmplifiedSearchRadiusExtension"].GetDouble()

        if 'MaxAmplificationRatioOfSearchRadius' in DEM_parameters.keys():
            self.max_amplification_ratio_of_search_radius = DEM_parameters["MaxAmplificationRatioOfSearchRadius"].GetDouble()
        else:
            self.max_amplification_ratio_of_search_radius = 0.0

        if not "PostPoissonRatio" in DEM_parameters.keys():
            self.poisson_ratio_option = 0
        else:
            self.poisson_ratio_option = DEM_parameters["PostPoissonRatio"].GetBool()

        if not "PoissonEffectOption" in DEM_parameters.keys():
            self.poisson_effect_option = False
        else:
            self.poisson_effect_option = DEM_parameters["PoissonEffectOption"].GetBool()

        if not "ShearStrainParallelToBondOption" in DEM_parameters.keys():
            self.shear_strain_parallel_to_bond_option = False
        else:
            self.shear_strain_parallel_to_bond_option = DEM_parameters["ShearStrainParallelToBondOption"].GetBool()

        if (self.poisson_effect_option or self.shear_strain_parallel_to_bond_option):
            self.compute_stress_tensor_option = 1


    def CreateCPlusPlusStrategy(self):

        self.SetVariablesAndOptions()

        # ADDITIONAL VARIABLES AND OPTIONS
        self.spheres_model_part.ProcessInfo.SetValue(AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION, self.amplified_continuum_search_radius_extension)
        self.spheres_model_part.ProcessInfo.SetValue(MAX_AMPLIFICATION_RATIO_OF_THE_SEARCH_RADIUS, self.max_amplification_ratio_of_search_radius)
        if self.contact_mesh_option:
            self.spheres_model_part.ProcessInfo.SetValue(CONTACT_MESH_OPTION, 1)
        else:
            self.spheres_model_part.ProcessInfo.SetValue(CONTACT_MESH_OPTION, 0)

        if ((self.test_type == "Triaxial") or (self.test_type == "Hydrostatic")):
            self.spheres_model_part.ProcessInfo.SetValue(TRIAXIAL_TEST_OPTION, 1)
        else:
            self.spheres_model_part.ProcessInfo.SetValue(TRIAXIAL_TEST_OPTION, 0)

        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, POISSON_EFFECT_OPTION, self.poisson_effect_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, SHEAR_STRAIN_PARALLEL_TO_BOND_OPTION, self.shear_strain_parallel_to_bond_option)

        for properties in self.spheres_model_part.Properties:
            ContinuumConstitutiveLawString = properties[DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME]
            ContinuumConstitutiveLaw = globals().get(ContinuumConstitutiveLawString)()
            if ContinuumConstitutiveLaw.CheckRequirementsOfStressTensor():
                self.spheres_model_part.ProcessInfo.SetValue(COMPUTE_STRESS_TENSOR_OPTION, 1)
                break

        strategy_parameters = self.DEM_parameters["strategy_parameters"]

        if (self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Velocity_Verlet'):
            self.cplusplus_strategy = ContinuumVelocityVerletSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                                            self.delta_option, self.creator_destructor, self.dem_fem_search, self.search_strategy, strategy_parameters)
        else:
            self.cplusplus_strategy = ContinuumExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                  self.delta_option, self.creator_destructor, self.dem_fem_search, self.search_strategy, strategy_parameters)

    def Initialize(self):
        self.cplusplus_strategy.Initialize()  # Calls the cplusplus_strategy Initialize function (initializes all elements and performs other necessary tasks before starting the time loop) (C++)

    def SetContinuumType(self):
        self.continuum_type = True

    def Initial_Critical_Time(self):        # Calls deprecated function
        (self.cplusplus_strategy).InitialTimeStepCalculation()

    def PrepareContactElementsForPrinting(self):
        (self.cplusplus_strategy).PrepareContactElementsForPrinting()

    def AddAdditionalVariables(self, spheres_model_part, DEM_parameters):
        spheres_model_part.AddNodalSolutionStepVariable(COHESIVE_GROUP)  # Continuum group
        spheres_model_part.AddNodalSolutionStepVariable(SKIN_SPHERE)

    def ModifyProperties(self, properties, param = 0):
        BaseExplicitStrategy.ModifyProperties(self, properties, param)

        if not param:
            ContinuumConstitutiveLawString = properties[DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME]
            ContinuumConstitutiveLaw = globals().get(ContinuumConstitutiveLawString)()
            ContinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties, True)
