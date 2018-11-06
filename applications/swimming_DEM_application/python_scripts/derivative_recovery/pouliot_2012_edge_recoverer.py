from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
from . import recoverer


class Pouliot2012EdgeDerivativesRecoverer(recoverer.DerivativesRecoverer):
    def GetFieldUtility(self):
        import math
        a = math.pi / 4
        d = math.pi / 2*0

        self.flow_field = EthierFlowField(a, d)
        space_time_set = SpaceTimeSet()
        self.field_utility = FluidFieldUtility(space_time_set, self.flow_field, 1000.0, 1e-6)
        return self.field_utility
    def __init__(self, pp, model_part):
        recoverer.DerivativesRecoverer.__init__(self, pp, model_part)
        self.dimension = pp.domain_size
        self.model_part = model_part
        self.use_lumped_mass_matrix = pp.CFD_DEM["material_acceleration_calculation_type"].GetInt() == 3
        self.recovery_model_part = ModelPart("PostGradientFluidPart")
        self.custom_functions_tool = CustomFunctionsCalculator3D()
        self.calculate_vorticity = pp.CFD_DEM["vorticity_calculation_type"].GetInt() > 0
        self.GetFieldUtility()
        self.CreateCPluPlusStrategies()

    def FillUpModelPart(self, element_type):
        for node in self.model_part.Nodes:
            self.recovery_model_part.AddNode(node, 0)
        self.recovery_model_part.ProcessInfo = self.model_part.ProcessInfo
        self.meshing_tool = self.GetMeshingTool()
        self.meshing_tool.FillUpEdgesModelPartFromSimplicesModelPart(self.recovery_model_part, self.model_part, element_type)

    def FillSetOfAllEdges(self, set_of_all_edges):
        for elem in self.model_part.Elements:
            for i, first_node in enumerate(elem.GetNodes()[: - 1]):
                for j, second_node in enumerate(elem.GetNodes()[i + 1 :]):
                    edge_ids = tuple(sorted((first_node.Id, second_node.Id)))
                    set_of_all_edges.add(edge_ids)

    def GetMeshingTool(self):
        if self.dimension == 2:
            return DerivativeRecoveryMeshingTools2D()
        else:
            return DerivativeRecoveryMeshingTools3D()

    def CreateCPluPlusStrategies(self, echo_level = 1):
        from KratosMultiphysics.ExternalSolversApplication import SuperLUIterativeSolver
        # from KratosMultiphysics.ExternalSolversApplication import SuperLUSolver
        # linear_solver = SuperLUIterativeSolver()
        scheme = ResidualBasedIncrementalUpdateStaticScheme()
        amgcl_smoother = AMGCLSmoother.SPAI0
        amgcl_krylov_type = AMGCLIterativeSolverType.BICGSTAB_WITH_GMRES_FALLBACK
        tolerance = 1e-8
        max_iterations = 400
        verbosity = 2 #0->shows no information, 1->some information, 2->all the information
        gmres_size = 400

        # if self.use_lumped_mass_matrix:
        #     linear_solver = CGSolver()
        # else:
        linear_solver = AMGCLSolver(amgcl_smoother, amgcl_krylov_type, tolerance, max_iterations, verbosity,gmres_size)
        # linear_solver = SuperLUIterativeSolver()
        # linear_solver = CGSolver()
        # linear_solver = SkylineLUFactorizationSolver()
        # linear_solver = SuperLUSolver()
        # linear_solver = ITSOL_ARMS_Solver()
        # linear_solver = MKLPardisoSolver()
        # linear_solver = AMGCLSolver(amgcl_smoother, amgcl_krylov_type, tolerance, max_iterations, verbosity,gmres_size)
        self.recovery_strategy = ResidualBasedDerivativeRecoveryStrategy(self.recovery_model_part, scheme, linear_solver, False, False, False, False)

        self.recovery_strategy.SetEchoLevel(echo_level)

    def AddDofs(self, DOF_variables):
        for node in self.recovery_model_part.Nodes:
            for var in DOF_variables:
                node.AddDof(var)

        print("DOFs for the derivative recovery solvers added correctly")

    def Solve(self):
        pass

    def SetToZero(self, variable):
        if type(variable).__name__ == 'DoubleVariable':
            self.custom_functions_tool.SetValueOfAllNotes(self.recovery_model_part, 0.0, variable)
        elif type(variable).__name__ == 'Array1DVariable3':
            self.custom_functions_tool.SetValueOfAllNotes(self.recovery_model_part, ZeroVector(3), variable)

class Pouliot2012EdgeGradientRecoverer(Pouliot2012EdgeDerivativesRecoverer, recoverer.VorticityRecoverer):
    def __init__(self, pp, model_part):
        Pouliot2012EdgeDerivativesRecoverer.__init__(self, pp, model_part)
        self.element_type = self.GetElementType(pp)
        self.FillUpModelPart(self.element_type)
        self.DOFs = self.GetDofs(pp)
        self.AddDofs(self.DOFs)
        self.calculate_vorticity = self.pp.CFD_DEM["lift_force_type"].GetInt()

    def Solve(self):
        print("\nSolving for the fluid acceleration...")
        sys.stdout.flush()
        self.SetToZero(VELOCITY_COMPONENT_GRADIENT)
        self.recovery_strategy.Solve()

    def RecoverGradientOfVelocity(self):
        self.model_part.ProcessInfo[CURRENT_COMPONENT] = 0
        self.Solve()
        self.custom_functions_tool.CopyValuesFromFirstToSecond(self.model_part, VELOCITY_COMPONENT_GRADIENT, VELOCITY_X_GRADIENT)
        self.model_part.ProcessInfo[CURRENT_COMPONENT] = 1
        self.Solve()
        self.custom_functions_tool.CopyValuesFromFirstToSecond(self.model_part, VELOCITY_COMPONENT_GRADIENT, VELOCITY_Y_GRADIENT)
        self.model_part.ProcessInfo[CURRENT_COMPONENT] = 2
        self.Solve() # and there is no need to copy anything
        self.custom_functions_tool.CopyValuesFromFirstToSecond(self.model_part, VELOCITY_COMPONENT_GRADIENT, VELOCITY_Z_GRADIENT)

        if self.calculate_vorticity:
            self.cplusplus_recovery_tool.CalculateVorticityFromGradient(self.model_part, VELOCITY_X_GRADIENT, VELOCITY_Y_GRADIENT, VELOCITY_Z_GRADIENT, VORTICITY)

    def RecoverGradientOfVelocityComponent(self, component):
        self.model_part.ProcessInfo[CURRENT_COMPONENT] = component
        self.Solve()
        if self.calculate_vorticity:
            self.cplusplus_recovery_tool.CalculateVorticityContributionOfTheGradientOfAComponent(self.model_part, VELOCITY_COMPONENT_GRADIENT, VORTICITY)

    def GetElementType(self, pp):
        if pp.domain_size == 2:
            return 'ComputeGradientPouliot20122DEdge'
        else:
            return 'ComputeGradientPouliot20123DEdge'

    def GetDofs(self, pp):
        if pp.domain_size == 2:
            return (VELOCITY_COMPONENT_GRADIENT_X, VELOCITY_COMPONENT_GRADIENT_Y)
        else:
            return (VELOCITY_COMPONENT_GRADIENT_X, VELOCITY_COMPONENT_GRADIENT_Y, VELOCITY_COMPONENT_GRADIENT_Z)

class Pouliot2012EdgeMaterialAccelerationRecoverer(Pouliot2012EdgeGradientRecoverer, recoverer.MaterialAccelerationRecoverer):
    def __init__(self, pp, model_part):
        Pouliot2012EdgeGradientRecoverer.__init__(self, pp, model_part)
        self.store_full_gradient = self.pp.CFD_DEM["store_full_gradient_option"].GetBool()

    def RecoverMaterialAcceleration(self):
        if self.store_full_gradient:
            self.RecoverGradientOfVelocity()
            self.RecoverMaterialAccelerationFromGradient()
        else:
            self.RecoverGradientOfVelocityComponent(0)
            self.cplusplus_recovery_tool.CalculateVectorMaterialDerivativeComponent(self.model_part, VELOCITY_COMPONENT_GRADIENT, ACCELERATION, MATERIAL_ACCELERATION)
            self.RecoverGradientOfVelocityComponent(1)
            self.cplusplus_recovery_tool.CalculateVectorMaterialDerivativeComponent(self.model_part, VELOCITY_COMPONENT_GRADIENT, ACCELERATION, MATERIAL_ACCELERATION)
            self.RecoverGradientOfVelocityComponent(2)
            self.cplusplus_recovery_tool.CalculateVectorMaterialDerivativeComponent(self.model_part, VELOCITY_COMPONENT_GRADIENT, ACCELERATION, MATERIAL_ACCELERATION)

class Pouliot2012EdgeLaplacianRecoverer(Pouliot2012EdgeMaterialAccelerationRecoverer, recoverer.LaplacianRecoverer):
    def __init__(self, pp, model_part):
        Pouliot2012EdgeDerivativesRecoverer.__init__(self, pp, model_part)
        self.element_type = self.GetElementType(pp)
        self.condition_type = self.GetConditionType(pp)
        self.FillUpModelPart(self.element_type)
        self.DOFs = self.GetDofs(pp)
        self.AddDofs(self.DOFs)

    def RecoverVelocityLaplacian(self):
        print("\nSolving for the laplacian...")
        self.SetToZero(VELOCITY_LAPLACIAN)
        self.recovery_strategy.Solve()

    def GetElementType(self, pp):
        if pp.domain_size == 2:
            return 'ComputeVelocityLaplacianSimplex2D'
        else:
            return 'ComputeVelocityLaplacianSimplex3D'

    def GetConditionType(self, pp):
        if pp.domain_size == 2:
            return 'ComputeLaplacianSimplexCondition2D'
        else:
            return 'ComputeLaplacianSimplexCondition3D'

    def GetDofs(self, pp):
        if pp.domain_size == 2:
            return (VELOCITY_LAPLACIAN_X, VELOCITY_LAPLACIAN_Y)
        else:
            return (VELOCITY_LAPLACIAN_X, VELOCITY_LAPLACIAN_Y, VELOCITY_LAPLACIAN_Z)
