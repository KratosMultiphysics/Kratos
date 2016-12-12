from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()
import sys

class DerivativeRecoveryStrategy:
    def __init__(self, pp, fluid_model_part, derivative_recovery_tool = None, custom_functions_tool = None):
        self.fluid_model_part = fluid_model_part
        self.derivative_recovery_tool = derivative_recovery_tool
        self.custom_functions_tool = custom_functions_tool
        self.pp = pp
        model_part_cloner = ConnectivityPreserveModeler()
        self.do_recover_acceleration = False
        self.do_recover_laplacian = False

        if pp.CFD_DEM.material_acceleration_calculation_type == 3:
            self.do_recover_acceleration = True
            self.acc_model_part = ModelPart("PostAccelerationFluidPart")
            model_part_cloner.GenerateModelPart(fluid_model_part, self.acc_model_part, "ComputeMaterialDerivativeSimplex3D", "ComputeLaplacianSimplexCondition3D")

        if pp.CFD_DEM.laplacian_calculation_type == 3:
            self.do_recover_laplacian = True
            self.lapl_model_part = ModelPart("PostLaplacianFluidPart")
            model_part_cloner.GenerateModelPart(fluid_model_part, self.lapl_model_part, "ComputeLaplacianSimplex3D", "ComputeLaplacianSimplexCondition3D")

        self.AddDofs(fluid_model_part)
        self.CreateCPluPlusStrategies(pp.CFD_DEM.recovery_echo_level)

    def AddDofs(self, model_part, config = None):
        if self.do_recover_acceleration:
            for node in self.acc_model_part.Nodes:
                node.AddDof(MATERIAL_ACCELERATION_X)
                node.AddDof(MATERIAL_ACCELERATION_Y)
                node.AddDof(MATERIAL_ACCELERATION_Z)
        if self.do_recover_laplacian:
            for node in self.lapl_model_part.Nodes:
                node.AddDof(VELOCITY_LAPLACIAN_X)
                node.AddDof(VELOCITY_LAPLACIAN_Y)
                node.AddDof(VELOCITY_LAPLACIAN_Z)

        print("dofs for the derivative recovery solvers added correctly")

    def CreateCPluPlusStrategies(self, echo_level = 1):
        linear_solver = CGSolver()
        scheme = ResidualBasedIncrementalUpdateStaticScheme()

        if self.do_recover_acceleration:
            self.acc_strategy = ResidualBasedLinearStrategy(self.acc_model_part, scheme, linear_solver, False, True, False, False)
            self.acc_strategy.SetEchoLevel(echo_level)

        if self.do_recover_laplacian:
            self.lapl_strategy = ResidualBasedLinearStrategy(self.lapl_model_part, scheme, linear_solver, False, True, False, False)
            self.lapl_strategy.SetEchoLevel(echo_level)

    def Solve(self):
        if self.do_recover_acceleration:
            print("\nSolving for the fluid acceleration...")
            sys.stdout.flush()

            for node in self.fluid_model_part.Nodes:
                node.SetSolutionStepValue(MATERIAL_ACCELERATION_X, 0.0)
                node.SetSolutionStepValue(MATERIAL_ACCELERATION_Y, 0.0)
                node.SetSolutionStepValue(MATERIAL_ACCELERATION_Z, 0.0)

            self.acc_strategy.Solve()

            print("\nFinished solving for the fluid acceleratio.\n")
            sys.stdout.flush()

        if self.do_recover_laplacian:
            print("\nSolving for the velocity laplacian...")
            sys.stdout.flush()

            for node in self.fluid_model_part.Nodes:
                node.SetSolutionStepValue(VELOCITY_LAPLACIAN_X, 0.0)
                node.SetSolutionStepValue(VELOCITY_LAPLACIAN_Y, 0.0)
                node.SetSolutionStepValue(VELOCITY_LAPLACIAN_Z, 0.0)

            self.lapl_strategy.Solve()

            print("\nFinished solving for the velocity laplacian.\n")
            sys.stdout.flush()

    def Recover(self):
        if self.pp.CFD_DEM.gradient_calculation_type == 1:
            self.custom_functions_tool.CalculatePressureGradient(self.fluid_model_part)
        elif self.pp.CFD_DEM.gradient_calculation_type == 2:
            self.derivative_recovery_tool.RecoverSuperconvergentGradient(self.fluid_model_part, PRESSURE, PRESSURE_GRADIENT)
        if self.pp.CFD_DEM.material_acceleration_calculation_type == 1:
            self.derivative_recovery_tool.CalculateVectorMaterialDerivative(self.fluid_model_part, VELOCITY, ACCELERATION, MATERIAL_ACCELERATION)
        elif self.pp.CFD_DEM.material_acceleration_calculation_type == 2 and self.pp.CFD_DEM.laplacian_calculation_type == 2:
            self.derivative_recovery_tool.RecoverSuperconvergentMatDerivAndLaplacian(fluid_model_part, VELOCITY, ACCELERATION, MATERIAL_ACCELERATION, VELOCITY_LAPLACIAN)
        else:
            if self.pp.CFD_DEM.laplacian_calculation_type == 1:
                self.derivative_recovery_tool.CalculateVectorLaplacian(self.fluid_model_part, VELOCITY, VELOCITY_LAPLACIAN)
            elif self.pp.CFD_DEM.laplacian_calculation_type == 2:
                self.derivative_recovery_tool.RecoverSuperconvergentLaplacian(self.fluid_model_part, VELOCITY, VELOCITY_LAPLACIAN)
            elif self.pp.CFD_DEM.material_acceleration_calculation_type == 2:
                self.derivative_recovery_tool.RecoverSuperconvergentMatDeriv(self.fluid_model_part, VELOCITY, ACCELERATION, MATERIAL_ACCELERATION)
        if self.pp.CFD_DEM.material_acceleration_calculation_type == 3 or self.pp.CFD_DEM.laplacian_calculation_type == 3:
            self.Solve()


