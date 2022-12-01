from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import KratosMultiphysics.DEMApplication.sphere_strategy as SolverStrategy
BaseExplicitStrategy = SolverStrategy.ExplicitStrategy

import math

class ExplicitStrategy(BaseExplicitStrategy):

    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures):

        BaseExplicitStrategy.__init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures)

        self.print_skin_sphere = 0 #TODO: check if this variable is important. There's a similar one in DEM_procedures called PostSkinSphere
        if "PostSkinSphere" in DEM_parameters.keys():
            self.print_skin_sphere = DEM_parameters["PostSkinSphere"].GetBool()

        if "DontSearchUntilFailure" in DEM_parameters.keys(): #TODO: important Todo. When Json gets divided in encapsulated parts, all these checks should be done in one functions, comparing with defaults!
            if DEM_parameters["DontSearchUntilFailure"].GetBool():
                print ("Search is not active until a bond is broken.")
                self.search_control = 0
                if (len(fem_model_part.Nodes) > 0 or DEM_parameters["TestType"].GetString() == "BTS"):   #MSI. This activates the search since there are fem contact elements. however only the particle - fem search should be active.
                    Logger.PrintWarning("DEM", "WARNING!: Search should be activated since there might contact with FEM.")

        if not "TestType" in DEM_parameters.keys():
            self.test_type = "None"
        else:
            self.test_type = DEM_parameters["TestType"].GetString()

        self.continuum_search_radius_amplification_factor = DEM_parameters["AmplifiedSearchRadiusExtension"].GetDouble()

        if 'MaxAmplificationRatioOfSearchRadius' in DEM_parameters.keys():
            self.max_amplification_ratio_of_search_radius = DEM_parameters["MaxAmplificationRatioOfSearchRadius"].GetDouble()
        else:
            self.max_amplification_ratio_of_search_radius = 0.0

        self.local_coordination_number_option = DEM_parameters["LocalCoordinationNumberOption"].GetBool()
        self.global_coordination_number_option = DEM_parameters["GlobalCoordinationNumberOption"].GetBool()

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

        if not "MaxNumberOfIntactBondsToConsiderASphereBroken" in DEM_parameters.keys():
            self.max_number_of_intact_bonds_to_consider_a_sphere_broken = 0
        else:
            self.max_number_of_intact_bonds_to_consider_a_sphere_broken = DEM_parameters["MaxNumberOfIntactBondsToConsiderASphereBroken"].GetDouble()

        if not "AutomaticSkinComputation" in DEM_parameters.keys():
            self.automatic_skin_computation = False
        else:
            self.automatic_skin_computation = DEM_parameters["AutomaticSkinComputation"].GetBool()

        if not "SkinFactorRadius" in DEM_parameters.keys():
            self.skin_factor_radius = 1.0
        else:
            self.skin_factor_radius = DEM_parameters["SkinFactorRadius"].GetDouble()

    def CreateCPlusPlusStrategy(self):

        self.SetVariablesAndOptions()

        # ADDITIONAL VARIABLES AND OPTIONS
        self.spheres_model_part.ProcessInfo.SetValue(CONTINUUM_SEARCH_RADIUS_AMPLIFICATION_FACTOR, self.continuum_search_radius_amplification_factor)
        self.spheres_model_part.ProcessInfo.SetValue(MAX_AMPLIFICATION_RATIO_OF_THE_SEARCH_RADIUS, self.max_amplification_ratio_of_search_radius)
        self.spheres_model_part.ProcessInfo.SetValue(LOCAL_COORDINATION_NUMBER_OPTION, self.local_coordination_number_option)
        self.spheres_model_part.ProcessInfo.SetValue(GLOBAL_COORDINATION_NUMBER_OPTION, self.global_coordination_number_option)

        if ((self.test_type == "Triaxial") or (self.test_type == "Hydrostatic")):
            self.spheres_model_part.ProcessInfo.SetValue(TRIAXIAL_TEST_OPTION, 1)
        else:
            self.spheres_model_part.ProcessInfo.SetValue(TRIAXIAL_TEST_OPTION, 0)

        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, POISSON_EFFECT_OPTION, self.poisson_effect_option)
        self.SetOneOrZeroInProcessInfoAccordingToBoolValue(self.spheres_model_part, SHEAR_STRAIN_PARALLEL_TO_BOND_OPTION, self.shear_strain_parallel_to_bond_option)
        self.spheres_model_part.ProcessInfo.SetValue(MAX_NUMBER_OF_INTACT_BONDS_TO_CONSIDER_A_SPHERE_BROKEN, self.max_number_of_intact_bonds_to_consider_a_sphere_broken)
        self.spheres_model_part.ProcessInfo.SetValue(AUTOMATIC_SKIN_COMPUTATION, self.automatic_skin_computation)
        self.spheres_model_part.ProcessInfo.SetValue(SKIN_FACTOR_RADIUS, self.skin_factor_radius)

        for properties in self.spheres_model_part.Properties:
            for subproperties in properties.GetSubProperties():
                if subproperties.Has(DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME):
                    continuum_constitutive_law_name = subproperties[DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME]
                    continuum_constitutive_law_instance = globals().get(continuum_constitutive_law_name)()
                    if continuum_constitutive_law_instance.CheckRequirementsOfStressTensor():
                        self.spheres_model_part.ProcessInfo.SetValue(COMPUTE_STRESS_TENSOR_OPTION, 1)
                        break

        if (self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Velocity_Verlet'):
            self.cplusplus_strategy = ContinuumVelocityVerletSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                                            self.delta_option, self.creator_destructor, self.dem_fem_search, self.search_strategy, self.solver_settings)
        else:
            self.cplusplus_strategy = ContinuumExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                  self.delta_option, self.creator_destructor, self.dem_fem_search, self.search_strategy, self.solver_settings)

    def BeforeInitialize(self):
        self.CreateCPlusPlusStrategy()
        self.RebuildListOfDiscontinuumSphericParticles()
        self.RebuildListOfContinuumSphericParticles()
        self.SetNormalRadiiOnAllParticles()
        self.SetSearchRadiiOnAllParticles()

    def Initialize(self):
        self.cplusplus_strategy.Initialize()  # Calls the cplusplus_strategy Initialize function (initializes all elements and performs other necessary tasks before starting the time loop) (C++)

    def SetContinuumType(self):
        self.continuum_type = True






    def AddAdditionalVariables(self, spheres_model_part, DEM_parameters):
        spheres_model_part.AddNodalSolutionStepVariable(COHESIVE_GROUP)  # Continuum group
        spheres_model_part.AddNodalSolutionStepVariable(SKIN_SPHERE)

    def RebuildListOfContinuumSphericParticles(self):
        self.cplusplus_strategy.RebuildListOfContinuumSphericParticles()



