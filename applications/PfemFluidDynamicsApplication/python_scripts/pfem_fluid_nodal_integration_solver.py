from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
import KratosMultiphysics.DelaunayMeshingApplication  as KratosDelaunay

from KratosMultiphysics.PfemFluidDynamicsApplication import pfem_fluid_solver as BaseSolver

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

        physics_type = self.settings["physics_type"].GetString()

        if( physics_type == "fsi" ):
            self.fluid_solver = KratosPfemFluid.NodalTwoStepVPStrategyForFSI(self.computing_model_part,
                                                                        self.velocity_linear_solver,
                                                                        self.pressure_linear_solver,
                                                                        self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                        self.settings["velocity_tolerance"].GetDouble(),
                                                                        self.settings["pressure_tolerance"].GetDouble(),
                                                                        self.settings["maximum_pressure_iterations"].GetInt(),
                                                                        self.settings["time_order"].GetInt(),
                                                                        self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION])
        else:
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
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PRESSURE_ACCELERATION)


        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_ERROR_XX)

        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_VOLUME)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_CAUCHY_STRESS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_DEVIATORIC_CAUCHY_STRESS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_SFD_NEIGHBOURS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_SFD_NEIGHBOURS_ORDER)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_DEFORMATION_GRAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_DEFORMATION_GRAD_VEL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_SPATIAL_DEF_RATE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_VOLUMETRIC_DEF_RATE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_EQUIVALENT_STRAIN_RATE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_MEAN_MESH_SIZE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_TAU)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.NODAL_FREESURFACE_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.VOLUMETRIC_COEFFICIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.DEVIATORIC_COEFFICIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.INTERFACE_NODE)

        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_NODAL_VOLUME)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_NODAL_CAUCHY_STRESS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_NODAL_SFD_NEIGHBOURS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_NODAL_SFD_NEIGHBOURS_ORDER)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_NODAL_DEFORMATION_GRAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_NODAL_DEFORMATION_GRAD_VEL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_NODAL_SPATIAL_DEF_RATE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_NODAL_VOLUMETRIC_DEF_RATE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_NODAL_EQUIVALENT_STRAIN_RATE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_NODAL_MEAN_MESH_SIZE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_NODAL_TAU)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_NODAL_FREESURFACE_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_VOLUMETRIC_COEFFICIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_DEVIATORIC_COEFFICIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.SOLID_INTERFACE_NODE)

        # Pfem Extra Vars
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONTACT_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONTACT_NORMAL)

        self.main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.OFFSET)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.SHRINK_FACTOR)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.MEAN_ERROR)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDelaunay.RIGID_WALL)

        self.main_model_part.AddNodalSolutionStepVariable(KratosPfemFluid.PROPERTY_ID)

        print("::[Pfem Fluid Solver]:: Variables ADDED")


    def InitializeSolutionStep(self):
        #self.fluid_solver.InitializeSolutionStep()
        if self._TimeBufferIsInitialized():
            self.fluid_solver.InitializeSolutionStep()

        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            adaptive_time_interval = KratosPfemFluid.AdaptiveTimeIntervalProcess(self.main_model_part,self.settings["echo_level"].GetInt())
            adaptive_time_interval.Execute()

        #pass
        #unactive_peak_elements = False
        #unactive_sliver_elements = False
        #set_active_flag = KratosPfemFluid.SetActiveFlagProcess(self.main_model_part,unactive_peak_elements,unactive_sliver_elements,self.settings["echo_level"].GetInt())
        #set_active_flag.Execute()

        # split_elements = KratosPfemFluid.SplitElementsProcess(self.main_model_part,self.settings["echo_level"].GetInt())
        # split_elements.ExecuteInitialize()

#


