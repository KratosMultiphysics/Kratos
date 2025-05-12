# importing the Kratos Library
import KratosMultiphysics as Kratos
from KratosMultiphysics import Vector, ModelPart
import KratosMultiphysics.SwimmingDEMApplication as SDEM
from . import recoverer
import KratosMultiphysics.SwimmingDEMApplication.parameters_tools as PT

import sys


class L2ProjectionDerivativesRecoverer(recoverer.DerivativesRecoverer):
    def __init__(self, project_parameters, model_part):
        recoverer.DerivativesRecoverer.__init__(self, project_parameters, model_part)

        self.model_part = model_part
        self.recovery_model_part = Kratos.Model().CreateModelPart("PostGradientFluidPart")
        self.custom_functions_tool = SDEM.CustomFunctionsCalculator3D()
        self.calculate_vorticity = (project_parameters["vorticity_calculation_type"].GetInt() > 0
                                    or PT.RecursiveFindParametersWithCondition(project_parameters["properties"],
                                                                               'vorticity_induced_lift_parameters'))

        if self.use_lumped_mass_matrix:
            self.model_part.ProcessInfo[Kratos.COMPUTE_LUMPED_MASS_MATRIX] = 1
        else:
            self.model_part.ProcessInfo[Kratos.COMPUTE_LUMPED_MASS_MATRIX] = 0
        self.CreateCPluPlusStrategies()

        self.num_points_element = model_part.GetElement(1).GetGeometry().PointsNumber()
        self.num_points_condition = 4 if self.num_points_element == 3 else 6

    def FillUpModelPart(self, element_type, condition_type):
        model_part_cloner = Kratos.ConnectivityPreserveModeler()
        model_part_cloner.GenerateModelPart(self.model_part, self.recovery_model_part, element_type, condition_type)

    def CreateCPluPlusStrategies(self, echo_level=1):
        try:
            from Kratos.ExternalSolversApplication import SuperLUIterativeSolver
            linear_solver = SuperLUIterativeSolver()
        except:
            linear_solver = Kratos.SkylineLUFactorizationSolver()
        # from KratosMultiphysics.ExternalSolversApplication import SuperLUIterativeSolver
        # linear_solver = SuperLUIterativeSolver()
        # from KratosMultiphysics.ExternalSolversApplication import SuperLUSolver
        scheme = Kratos.ResidualBasedIncrementalUpdateStaticScheme()
        amgcl_smoother = Kratos.AMGCLSmoother.ILU0
        amgcl_krylov_type = Kratos.AMGCLIterativeSolverType.BICGSTAB2
        tolerance = 1e-12
        max_iterations = 1000
        verbosity = 0  # 0->shows no information, 1->some information, 2->all the information
        gmres_size = 50

        if not self.use_lumped_mass_matrix:
            linear_solver = Kratos.AMGCLSolver(amgcl_smoother, amgcl_krylov_type, tolerance, max_iterations, verbosity, gmres_size)

        builder_and_solver = Kratos.ResidualBasedBlockBuilderAndSolver(linear_solver)
        self.recovery_strategy = SDEM.ResidualBasedDerivativeRecoveryStrategy(self.recovery_model_part, scheme, linear_solver, builder_and_solver, False, True, False, False)
        self.recovery_strategy.SetEchoLevel(echo_level)

    def AddDofs(self, DOF_variables):
        if len(self.recovery_model_part.Nodes) == 0:
            print(f"Warning: The derivative recovery model part has no nodes:\n{self.recovery_model_part}")
        for node in self.recovery_model_part.Nodes:
            for var in DOF_variables:
                node.AddDof(var)
        Kratos.Logger.PrintInfo("SwimmingDEM", "DOFs for the derivative recovery solvers added correctly")

    def Solve(self):
        pass

    def SetToZero(self, variable):
        if type(variable).__name__ == 'DoubleVariable':
            self.custom_functions_tool.SetValueOfAllNotes(self.model_part, 0.0, variable)
        elif type(variable).__name__ == 'Array1DVariable3':
            self.custom_functions_tool.SetValueOfAllNotes(self.model_part, Vector([0, 0, 0]), variable)


class L2ProjectionGradientRecoverer(L2ProjectionDerivativesRecoverer, recoverer.VorticityRecoverer):
    def __init__(self, project_parameters, model_part):
        self.use_lumped_mass_matrix = project_parameters["material_acceleration_calculation_type"].GetInt() == 3
        L2ProjectionDerivativesRecoverer.__init__(self, project_parameters, model_part)

        self.element_type = f"ComputeComponentGradientSimplex3D{self.num_points_element}N"
        self.condition_type = f"ComputeLaplacianSimplexCondition3D{self.num_points_condition}N"
        self.FillUpModelPart(self.element_type, self.condition_type)
        self.DOFs = (Kratos.VELOCITY_COMPONENT_GRADIENT_X, Kratos.VELOCITY_COMPONENT_GRADIENT_Y, Kratos.VELOCITY_COMPONENT_GRADIENT_Z)
        self.AddDofs(self.DOFs)
        self.calculate_vorticity = (project_parameters["vorticity_calculation_type"].GetInt() > 0
                                    or PT.RecursiveFindParametersWithCondition(project_parameters["properties"],
                                                                               'vorticity_induced_lift_parameters'))

    def Solve(self):
        Kratos.Logger.PrintInfo("SwimmingDEM", "(L2ProjectionGradientRecoverer)", "\nSolving for the fluid acceleration (L2)...")
        sys.stdout.flush()
        self.SetToZero(Kratos.VELOCITY_COMPONENT_GRADIENT)
        self.recovery_strategy.Solve()

    def RecoverGradientOfVelocity(self):
        self.model_part.ProcessInfo[Kratos.CURRENT_COMPONENT] = 0
        self.Solve()
        self.custom_functions_tool.CopyValuesFromFirstToSecond(self.model_part, Kratos.VELOCITY_COMPONENT_GRADIENT, Kratos.VELOCITY_X_GRADIENT)
        self.model_part.ProcessInfo[Kratos.CURRENT_COMPONENT] = 1
        self.Solve()
        self.custom_functions_tool.CopyValuesFromFirstToSecond(self.model_part, Kratos.VELOCITY_COMPONENT_GRADIENT, Kratos.VELOCITY_Y_GRADIENT)
        self.model_part.ProcessInfo[Kratos.CURRENT_COMPONENT] = 2
        self.Solve()
        self.custom_functions_tool.CopyValuesFromFirstToSecond(self.model_part, Kratos.VELOCITY_COMPONENT_GRADIENT, Kratos.VELOCITY_Z_GRADIENT)

        if self.calculate_vorticity:
            self.cplusplus_recovery_tool.CalculateVorticityFromGradient(self.model_part, Kratos.VELOCITY_X_GRADIENT, Kratos.VELOCITY_Y_GRADIENT, Kratos.VELOCITY_Z_GRADIENT, Kratos.VORTICITY)

    def RecoverGradientOfVelocityComponent(self, component):
        self.model_part.ProcessInfo[Kratos.CURRENT_COMPONENT] = component
        self.Solve()
        if self.calculate_vorticity:
            self.cplusplus_recovery_tool.CalculateVorticityContributionOfTheGradientOfAComponent(self.model_part, Kratos.VELOCITY_COMPONENT_GRADIENT, Kratos.VORTICITY)


class L2ProjectionMaterialAccelerationRecoverer(L2ProjectionGradientRecoverer, recoverer.MaterialAccelerationRecoverer):
    def __init__(self, project_parameters, model_part):
        self.use_lumped_mass_matrix = project_parameters["material_acceleration_calculation_type"].GetInt() == 3
        L2ProjectionGradientRecoverer.__init__(self, project_parameters, model_part)
        self.store_full_gradient = project_parameters["store_full_gradient_option"].GetBool()
        print("In constructor of L2ProjectionMaterialAccelerationRecoverer")

    def RecoverMaterialAcceleration(self):
        if self.store_full_gradient:
            # self.RecoverGradientOfVelocity()
            self.RecoverMaterialAccelerationFromGradient()
        else:
            self.RecoverGradientOfVelocityComponent(0)
            self.cplusplus_recovery_tool.CalculateVectorMaterialDerivativeComponent(self.model_part, Kratos.VELOCITY_COMPONENT_GRADIENT, Kratos.ACCELERATION, Kratos.MATERIAL_ACCELERATION)
            self.RecoverGradientOfVelocityComponent(1)
            self.cplusplus_recovery_tool.CalculateVectorMaterialDerivativeComponent(self.model_part, Kratos.VELOCITY_COMPONENT_GRADIENT, Kratos.ACCELERATION, Kratos.MATERIAL_ACCELERATION)
            self.RecoverGradientOfVelocityComponent(2)
            self.cplusplus_recovery_tool.CalculateVectorMaterialDerivativeComponent(self.model_part, Kratos.VELOCITY_COMPONENT_GRADIENT, Kratos.ACCELERATION, Kratos.MATERIAL_ACCELERATION)


class L2ProjectionDirectMaterialAccelerationRecoverer(L2ProjectionMaterialAccelerationRecoverer):
    def __init__(self, project_parameters, model_part):
        self.use_lumped_mass_matrix = project_parameters["material_acceleration_calculation_type"].GetInt() == 3
        L2ProjectionDerivativesRecoverer.__init__(self, project_parameters, model_part)
        self.element_type = f"ComputeMaterialDerivativeSimplex3D{self.num_points_element}N"
        self.condition_type = f"ComputeLaplacianSimplexCondition3D{self.num_points_condition}N"
        self.FillUpModelPart(self.element_type, self.condition_type)
        self.DOFs = (Kratos.MATERIAL_ACCELERATION_X, Kratos.MATERIAL_ACCELERATION_Y, Kratos.MATERIAL_ACCELERATION_Z)
        self.AddDofs(self.DOFs)

    def RecoverMaterialAcceleration(self):
        self.SetToZero(Kratos.MATERIAL_ACCELERATION)
        self.recovery_strategy.Solve()


class L2ProjectionFluidFractionGradientRecoverer(L2ProjectionDerivativesRecoverer):
    def __init__(self, project_parameters, model_part):
        self.use_lumped_mass_matrix = 0
        L2ProjectionDerivativesRecoverer.__init__(self, project_parameters, model_part)
        domain_size = model_part.ProcessInfo[Kratos.DOMAIN_SIZE]
        self.element_type = "ComputeFluidFractionGradient" + str(int(domain_size)) + "D" + str(int(self.num_points_element)) + "N"
        self.condition_type = "MonolithicWallCondition" + str(int(domain_size)) + "D" + str(int(self.num_points_condition)) + "N"
        self.FillUpModelPart(self.element_type, self.condition_type)
        self.DOFs = (Kratos.FLUID_FRACTION_GRADIENT_X, Kratos.FLUID_FRACTION_GRADIENT_Y, Kratos.FLUID_FRACTION_GRADIENT_Z)
        self.AddDofs(self.DOFs)

    def RecoverFluidFractionGradient(self, model_part):
        self.model_part = model_part
        self.element_type = "ComputeFluidFractionGradient2D4N"
        self.condition_type = "MonolithicWallCondition2D2N"
        self.recovery_model_part.ProcessInfo[Kratos.COMPUTE_LUMPED_MASS_MATRIX] == 1
        self.model_part.ProcessInfo[Kratos.COMPUTE_LUMPED_MASS_MATRIX] == 1
        self.FillUpModelPart(self.element_type, self.condition_type)
        self.DOFs = (Kratos.FLUID_FRACTION_GRADIENT_X, Kratos.FLUID_FRACTION_GRADIENT_Y, Kratos.FLUID_FRACTION_GRADIENT_Z)
        self.AddDofs(self.DOFs)
        self.CreateCPluPlusStrategies()
        self.SetToZero(Kratos.FLUID_FRACTION_GRADIENT)
        self.recovery_strategy.Solve()


class L2ProjectionLaplacianRecoverer(L2ProjectionMaterialAccelerationRecoverer, recoverer.LaplacianRecoverer):
    def __init__(self, project_parameters, model_part):
        self.use_lumped_mass_matrix = project_parameters["material_acceleration_calculation_type"].GetInt() == 3
        L2ProjectionDerivativesRecoverer.__init__(self, project_parameters, model_part)
        self.element_type = f"ComputeVelocityLaplacianSimplex3D{self.num_points_element}N"
        self.condition_type = f"ComputeLaplacianSimplexCondition3D{self.num_points_condition}N"
        self.FillUpModelPart(self.element_type, self.condition_type)
        self.DOFs = (Kratos.VELOCITY_LAPLACIAN_X, Kratos.VELOCITY_LAPLACIAN_Y, Kratos.VELOCITY_LAPLACIAN_Z)
        self.AddDofs(self.DOFs)

    def RecoverVelocityLaplacian(self):
        Kratos.Logger.PrintInfo("SwimmingDEM", "\nSolving for the laplacian...")
        self.SetToZero(Kratos.VELOCITY_LAPLACIAN)
        self.recovery_strategy.Solve()
