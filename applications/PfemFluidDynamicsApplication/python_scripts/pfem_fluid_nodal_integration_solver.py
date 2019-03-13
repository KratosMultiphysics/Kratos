from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

import pfem_fluid_solver as BaseSolver

def CreateSolver(model, parameters):
    return PfemFluidNodalIntegrationSolver(model, parameters)

class PfemFluidNodalIntegrationSolver(BaseSolver.PfemFluidSolver):

    def __init__(self, model, parameters):

        super(PfemFluidNodalIntegrationSolver, self).__init__(model, parameters)

        print("Construction of 2-step Pfem Fluid Nodal Integration Solver finished.")



    def Initialize(self):

        print("::[Pfem Fluid Nodal Integration Solver]:: -START-")

        # Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()

        self.fluid_solver = KratosPfemFluid.NodalTwoStepVPStrategy(self.computing_model_part,
                                                                   self.velocity_linear_solver,
                                                                   self.pressure_linear_solver,
                                                                   self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                   self.settings["velocity_tolerance"].GetDouble(),
                                                                   self.settings["pressure_tolerance"].GetDouble(),
                                                                   self.settings["maximum_pressure_iterations"].GetInt(),
                                                                   self.settings["time_order"].GetInt(),
                                                                   self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION])

        echo_level = self.settings["echo_level"].GetInt()

        # Set echo_level
        self.fluid_solver.SetEchoLevel(echo_level)

        # Set initialize flag
        if( self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True ):
            self.mechanical_solver.SetInitializePerformedFlag(True)


        # Check if everything is assigned correctly
        self.fluid_solver.Check()


        print("::[Pfem Fluid Nodal Integration Solver]:: -END- ")



    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BULK_MODULUS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POISSON_RATIO)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.YOUNG_MODULUS)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_ERROR)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)


        #VARIABLES FOR PAPANASTASIOU MODEL
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.FLOW_INDEX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.YIELD_SHEAR)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.ADAPTIVE_EXPONENT)

        #VARIABLES FOR MU-I RHEOLOGY MODEL
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.STATIC_FRICTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.DYNAMIC_FRICTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INERTIAL_NUMBER_ZERO)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DIAMETER)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.GRAIN_DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.REGULARIZATION_COEFFICIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INFINITE_FRICTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INERTIAL_NUMBER_ONE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.ALPHA_PARAMETER)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLUID_FRACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLUID_FRACTION_OLD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLUID_FRACTION_RATE)

        # PFEM fluid variables
        # self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NORMVELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.YIELDED)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.FREESURFACE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_ACCELERATION)

        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_VOLUME)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_CAUCHY_STRESS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_DEVIATORIC_CAUCHY_STRESS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_SFD_NEIGHBOURS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_SFD_NEIGHBOURS_ORDER)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_DEFORMATION_GRAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_DEFORMATION_GRAD_VEL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_SPATIAL_DEF_RATE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_SPATIAL_DEF_RATE_BIS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_VOLUMETRIC_DEF_RATE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_MEAN_MESH_SIZE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_FREESURFACE_AREA)

        print("::[Pfem Fluid Solver]:: Variables ADDED")


    def InitializeSolutionStep(self):
        #self.fluid_solver.InitializeSolutionStep()
        if self._TimeBufferIsInitialized():
            self.fluid_solver.InitializeSolutionStep()

        adaptive_time_interval = KratosPfemFluid.AdaptiveTimeIntervalProcess(self.main_model_part,self.settings["echo_level"].GetInt())
        adaptive_time_interval.Execute()

        #pass
        #unactive_peak_elements = False
        #unactive_sliver_elements = False
        #set_active_flag = KratosPfemFluid.SetActiveFlagProcess(self.main_model_part,unactive_peak_elements,unactive_sliver_elements,self.settings["echo_level"].GetInt())
        #set_active_flag.Execute()

        #split_elements = KratosPfemFluid.SplitElementsProcess(self.main_model_part,self.settings["echo_level"].GetInt())
        #split_elements.ExecuteInitialize()

#   Extra methods:: custom AFranci...
#
    def SetProperties(self):
        for el in self.main_model_part.Elements:
            density = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            viscosity = el.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
            bulk_modulus = el.Properties.GetValue(KratosMultiphysics.BULK_MODULUS)
            young_modulus = el.Properties.GetValue(KratosMultiphysics.YOUNG_MODULUS)
            poisson_ratio = el.Properties.GetValue(KratosMultiphysics.POISSON_RATIO)
            flow_index = el.Properties.GetValue(KratosPfemFluid.FLOW_INDEX)
            yield_shear = el.Properties.GetValue(KratosPfemFluid.YIELD_SHEAR)
            adaptive_exponent = el.Properties.GetValue(KratosPfemFluid.ADAPTIVE_EXPONENT)
            static_friction = elem.Properties.GetValue(KratosPfemFluid.STATIC_FRICTION)
            dynamic_friction = elem.Properties.GetValue(KratosPfemFluid.DYNAMIC_FRICTION)
            inertial_number_zero = elem.Properties.GetValue(KratosPfemFluid.INERTIAL_NUMBER_ZERO)
            grain_diameter = elem.Properties.GetValue(KratosPfemFluid.GRAIN_DIAMETER)
            grain_density = elem.Properties.GetValue(KratosPfemFluid.GRAIN_DENSITY)
            regularization_coefficient = elem.Properties.GetValue(KratosPfemFluid.REGULARIZATION_COEFFICIENT)
            inertial_number_one = elem.Properties.GetValue(KratosPfemFluid.INERTIAL_NUMBER_ONE)
            infinite_friction = elem.Properties.GetValue(KratosPfemFluid.INFINITE_FRICTION)
            alpha_parameter = elem.Properties.GetValue(KratosPfemFluid.ALPHA_PARAMETER)
            break

        print ("density: ",density)
        print ("viscosity: ",viscosity)
        print ("bulk_modulus: ",bulk_modulus)
        print ("young_modulus: ",young_modulus)
        print ("poisson_ratio: ",poisson_ratio)
        print ("flow_index: ",flow_index)
        print ("yield_shear: ",yield_shear)
        print ("adaptive_exponent: ",adaptive_exponent)
        print ("static_friction: ",static_friction)
        print ("dynamic_friction: ",dynamic_friction)
        print ("inertial_number_zero: ",inertial_number_zero)
        print ("grain_diameter: ",grain_diameter)
        print ("grain_density: ",grain_density)
        print ("regularization_coefficient: ",regularization_coefficient)
        print ("inertial_number_one: ",inertial_number_one)
        print ("infinite_friction: ",infinite_friction)
        print ("alpha_parameter: ",alpha_parameter)

#


#


