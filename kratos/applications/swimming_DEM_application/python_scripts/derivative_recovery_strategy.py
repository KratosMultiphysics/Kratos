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
        self.do_recover_gradient = False        
        self.do_recover_laplacian = False

        if pp.CFD_DEM.material_acceleration_calculation_type == 3 or pp.CFD_DEM.material_acceleration_calculation_type == 4:
            self.do_recover_acceleration = True
            self.acc_model_part = ModelPart("PostAccelerationFluidPart")
            model_part_cloner.GenerateModelPart(fluid_model_part, self.acc_model_part, "ComputeMaterialDerivativeSimplex3D", "ComputeLaplacianSimplexCondition3D")
        elif pp.CFD_DEM.material_acceleration_calculation_type == 5:
            self.do_recover_gradient = True
            self.acc_model_part = ModelPart("PostAccelerationFluidPart")
            model_part_cloner.GenerateModelPart(fluid_model_part, self.acc_model_part, "ComputeComponentGradientSimplex3D", "ComputeLaplacianSimplexCondition3D")

        if pp.CFD_DEM.laplacian_calculation_type == 3 or pp.CFD_DEM.laplacian_calculation_type == 4 or pp.CFD_DEM.laplacian_calculation_type == 5:
            self.do_recover_laplacian = True
            self.lapl_model_part = ModelPart("PostLaplacianFluidPart")
            model_part_cloner.GenerateModelPart(fluid_model_part, self.lapl_model_part, "ComputeLaplacianSimplex3D", "ComputeLaplacianSimplexCondition3D")

        self.AddDofs(fluid_model_part)
        self.CreateCPluPlusStrategies(pp.CFD_DEM.recovery_echo_level)
        
        if pp.CFD_DEM.material_acceleration_calculation_type == 3:
            fluid_model_part.ProcessInfo[COMPUTE_LUMPED_MASS_MATRIX] = 1
        elif pp.CFD_DEM.material_acceleration_calculation_type == 4 or pp.CFD_DEM.material_acceleration_calculation_type == 5:
            fluid_model_part.ProcessInfo[COMPUTE_LUMPED_MASS_MATRIX] = 0


    def AddDofs(self, model_part, config = None):
        if self.do_recover_acceleration:
            for node in self.acc_model_part.Nodes:
                node.AddDof(MATERIAL_ACCELERATION_X)
                node.AddDof(MATERIAL_ACCELERATION_Y)
                node.AddDof(MATERIAL_ACCELERATION_Z)
        if self.do_recover_gradient:
            for node in self.acc_model_part.Nodes:
                node.AddDof(VELOCITY_COMPONENT_GRADIENT_X)
                node.AddDof(VELOCITY_COMPONENT_GRADIENT_Y)
                node.AddDof(VELOCITY_COMPONENT_GRADIENT_Z)    

        if self.do_recover_laplacian:
            for node in self.lapl_model_part.Nodes:
                node.AddDof(VELOCITY_LAPLACIAN_X)
                node.AddDof(VELOCITY_LAPLACIAN_Y)
                node.AddDof(VELOCITY_LAPLACIAN_Z)

        print("dofs for the derivative recovery solvers added correctly")

    def CreateCPluPlusStrategies(self, echo_level = 1):
        from KratosMultiphysics.ExternalSolversApplication import SuperLUIterativeSolver
        #linear_solver = CGSolver()
        scheme = ResidualBasedIncrementalUpdateStaticScheme()        
        amgcl_smoother = AMGCLSmoother.ILU0
        amgcl_krylov_type = AMGCLIterativeSolverType.BICGSTAB
        tolerance = 1e-6
        max_iterations = 1000
        verbosity = 0 #0->shows no information, 1->some information, 2->all the information
        gmres_size = 50
        linear_solver = AMGCLSolver(amgcl_smoother,amgcl_krylov_type,tolerance,max_iterations,verbosity,gmres_size)

        if self.do_recover_acceleration or self.do_recover_gradient:
            self.acc_strategy = ResidualBasedLinearStrategy(self.acc_model_part, scheme, linear_solver, False, True, False, False)
            self.acc_strategy.SetEchoLevel(echo_level)

        if self.do_recover_laplacian:
            self.lapl_strategy = ResidualBasedLinearStrategy(self.lapl_model_part, scheme, linear_solver, False, True, False, False)
            self.lapl_strategy.SetEchoLevel(echo_level)

    def Solve(self, component = None):
        if self.do_recover_acceleration or self.do_recover_gradient:
            if self.do_recover_acceleration:
                print("\nSolving for the fluid acceleration...")
                sys.stdout.flush()

            if component == None:
                for node in self.fluid_model_part.Nodes:
                    node.SetSolutionStepValue(MATERIAL_ACCELERATION_X, 0.0)
                    node.SetSolutionStepValue(MATERIAL_ACCELERATION_Y, 0.0)
                    node.SetSolutionStepValue(MATERIAL_ACCELERATION_Z, 0.0)
            else:
                for node in self.fluid_model_part.Nodes:
                    if component == 0:
                        node.SetSolutionStepValue(MATERIAL_ACCELERATION_X, 0.0)
                    elif component == 1:
                        node.SetSolutionStepValue(MATERIAL_ACCELERATION_Y, 0.0)
                    elif component == 2:
                        node.SetSolutionStepValue(MATERIAL_ACCELERATION_Z, 0.0)   
                        
                    node.SetSolutionStepValue(VELOCITY_COMPONENT_GRADIENT_X, 0.0)
                    node.SetSolutionStepValue(VELOCITY_COMPONENT_GRADIENT_Y, 0.0)
                    node.SetSolutionStepValue(VELOCITY_COMPONENT_GRADIENT_Z, 0.0)                

            self.acc_strategy.Solve()
            
            if self.do_recover_acceleration:
                print("\nFinished solving for the fluid acceleration.\n")
                sys.stdout.flush()

        if component == None:                
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
        if self.pp.CFD_DEM.material_acceleration_calculation_type == 2 and self.pp.CFD_DEM.laplacian_calculation_type == 2:
            self.derivative_recovery_tool.RecoverSuperconvergentMatDerivAndLaplacian(self.fluid_model_part, VELOCITY, ACCELERATION, MATERIAL_ACCELERATION, VELOCITY_LAPLACIAN)
        else:
            if self.pp.CFD_DEM.laplacian_calculation_type == 1:
                self.derivative_recovery_tool.CalculateVectorLaplacian(self.fluid_model_part, VELOCITY, VELOCITY_LAPLACIAN)
            elif self.pp.CFD_DEM.laplacian_calculation_type == 2:
                self.derivative_recovery_tool.RecoverSuperconvergentLaplacian(self.fluid_model_part, VELOCITY, VELOCITY_LAPLACIAN)
            if self.pp.CFD_DEM.material_acceleration_calculation_type == 2:
                self.derivative_recovery_tool.RecoverSuperconvergentMatDeriv(self.fluid_model_part, VELOCITY, ACCELERATION, MATERIAL_ACCELERATION)
            elif self.pp.CFD_DEM.material_acceleration_calculation_type == 1:
                self.derivative_recovery_tool.CalculateVectorMaterialDerivative(self.fluid_model_part, VELOCITY, ACCELERATION, MATERIAL_ACCELERATION)
            elif self.pp.CFD_DEM.material_acceleration_calculation_type == 3 or self.pp.CFD_DEM.laplacian_calculation_type == 3 or self.pp.CFD_DEM.material_acceleration_calculation_type == 4 or self.pp.CFD_DEM.laplacian_calculation_type == 4:
                self.Solve()
            if self.pp.CFD_DEM.material_acceleration_calculation_type == 5:
                print("\nSolving for the fluid acceleration...")
                sys.stdout.flush()
                self.fluid_model_part.ProcessInfo[CURRENT_COMPONENT] = 0
                self.Solve(0)
                self.derivative_recovery_tool.CalculateVectorMaterialDerivativeComponent(self.fluid_model_part, VELOCITY_COMPONENT_GRADIENT, ACCELERATION, MATERIAL_ACCELERATION)         
                if self.pp.CFD_DEM.lift_force_type:
                    self.derivative_recovery_tool.CalculateVorticityContributionOfTheGradientOfAComponent(self.fluid_model_part, VELOCITY_COMPONENT_GRADIENT, VORTICITY)         
                self.fluid_model_part.ProcessInfo[CURRENT_COMPONENT] = 1
                self.Solve(1)
                self.derivative_recovery_tool.CalculateVectorMaterialDerivativeComponent(self.fluid_model_part, VELOCITY_COMPONENT_GRADIENT, ACCELERATION, MATERIAL_ACCELERATION)                
                if self.pp.CFD_DEM.lift_force_type:
                        self.derivative_recovery_tool.CalculateVorticityContributionOfTheGradientOfAComponent(self.fluid_model_part, VELOCITY_COMPONENT_GRADIENT, VORTICITY)                  
                self.fluid_model_part.ProcessInfo[CURRENT_COMPONENT] = 2
                self.Solve(2)
                self.derivative_recovery_tool.CalculateVectorMaterialDerivativeComponent(self.fluid_model_part, VELOCITY_COMPONENT_GRADIENT, ACCELERATION, MATERIAL_ACCELERATION)                
                if self.pp.CFD_DEM.lift_force_type:
                        self.derivative_recovery_tool.CalculateVorticityContributionOfTheGradientOfAComponent(self.fluid_model_part, VELOCITY_COMPONENT_GRADIENT, VORTICITY)                       
                print("\nFinished solving for the fluid acceleration.\n")
                sys.stdout.flush()


