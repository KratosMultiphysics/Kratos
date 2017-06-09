from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
from . import recoverer


class Pouliot2012EdgeDerivativesRecoverer(recoverer.DerivativesRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        recoverer.DerivativesRecoverer.__init__(self, pp, model_part, cplusplus_recovery_tool)
        self.model_part = model_part
        self.use_lumped_mass_matrix = pp.CFD_DEM.material_acceleration_calculation_type == 3
        self.recovery_model_part = ModelPart("PostGradientFluidPart")
        self.custom_functions_tool = CustomFunctionsCalculator3D()
        self.calculate_vorticity = pp.CFD_DEM.vorticity_calculation_type > 0

        self.CreateCPluPlusStrategies()

    def FillUpModelPart(self, element_type):
        set_of_all_edges = set()
        #self.FillSetOfAllEdges(set_of_all_edges)
        self.recovery_model_part.Nodes = self.model_part.Nodes
        self.recovery_model_part.ProcessInfo = self.model_part.ProcessInfo
        self.meshing_tool = DerivativeRecoveryMeshingTools()
        self.meshing_tool.FillUpEdgesModelPartFromTetrahedraModelPart(self.recovery_model_part, self.model_part, element_type)
        # for i, edge in enumerate(set_of_all_edges):
        #     self.recovery_model_part.CreateNewElement(element_type, i + 1000000, list(edge), self.model_part.GetProperties()[0])

    def FillSetOfAllEdges(self, set_of_all_edges):
        for elem in self.model_part.Elements:
            for i, first_node in enumerate(elem.GetNodes()[: - 1]):
                for j, second_node in enumerate(elem.GetNodes()[i + 1 :]):
                    edge_ids = tuple(sorted((first_node.Id, second_node.Id)))
                    set_of_all_edges.add(edge_ids)

    def CreateCPluPlusStrategies(self, echo_level = 1):
        #from KratosMultiphysics.ExternalSolversApplication import SuperLUIterativeSolver
        #linear_solver = SuperLUIterativeSolver()
        scheme = ResidualBasedIncrementalUpdateStaticScheme()
        amgcl_smoother = AMGCLSmoother.SPAI0
        amgcl_krylov_type = AMGCLIterativeSolverType.BICGSTAB_WITH_GMRES_FALLBACK
        tolerance = 1e-11
        max_iterations = 200
        verbosity = 2 #0->shows no information, 1->some information, 2->all the information
        gmres_size = 400

        if self.use_lumped_mass_matrix:
            linear_solver = CGSolver()
        else:
            linear_solver = AMGCLSolver(amgcl_smoother, amgcl_krylov_type, tolerance, max_iterations, verbosity,gmres_size)
        #linear_solver = SuperLUIterativeSolver()
        self.recovery_strategy = ResidualBasedDerivativeRecoveryStrategy(self.recovery_model_part, scheme, linear_solver, False, True, False, False)
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
            self.custom_functions_tool.SetValueOfAllNotes(self.model_part, 0.0, variable)
        elif type(variable).__name__ == 'Array1DVariable3':
            self.custom_functions_tool.SetValueOfAllNotes(self.model_part, ZeroVector(3), variable)

class Pouliot2012EdgeGradientRecoverer(Pouliot2012EdgeDerivativesRecoverer, recoverer.VorticityRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        Pouliot2012EdgeDerivativesRecoverer.__init__(self, pp, model_part, cplusplus_recovery_tool)
        self.element_type = "ComputeGradientPouliot20123DEdge"
        # self.condition_type = "Condition1D"
        self.FillUpModelPart(self.element_type)
        self.DOFs = (VELOCITY_Z_GRADIENT_X, VELOCITY_Z_GRADIENT_Y, VELOCITY_Z_GRADIENT_Z)
        self.AddDofs(self.DOFs)
        self.calculate_vorticity = self.pp.CFD_DEM.lift_force_type

    def Solve(self):
        print("\nSolving for the fluid acceleration...")
        sys.stdout.flush()
        self.SetToZero(VELOCITY_Z_GRADIENT)
        self.recovery_strategy.Solve()

    def RecoverGradientOfVelocity(self):
        self.model_part.ProcessInfo[CURRENT_COMPONENT] = 0
        self.Solve()
        self.custom_functions_tool.CopyValuesFromFirstToSecond(self.model_part, VELOCITY_Z_GRADIENT, VELOCITY_X_GRADIENT)
        self.model_part.ProcessInfo[CURRENT_COMPONENT] = 1
        self.Solve()
        self.custom_functions_tool.CopyValuesFromFirstToSecond(self.model_part, VELOCITY_Z_GRADIENT, VELOCITY_Y_GRADIENT)
        self.model_part.ProcessInfo[CURRENT_COMPONENT] = 2
        self.Solve() # and there is no need to copy anything
        if self.calculate_vorticity:
            self.cplusplus_recovery_tool.CalculateVorticityFromGradient(self.model_part, VELOCITY_X_GRADIENT, VELOCITY_Y_GRADIENT, VELOCITY_Z_GRADIENT, VORTICITY)

    def RecoverGradientOfVelocityComponent(self, component):
        self.model_part.ProcessInfo[CURRENT_COMPONENT] = component
        self.Solve()
        if self.calculate_vorticity:
            self.cplusplus_recovery_tool.CalculateVorticityContributionOfTheGradientOfAComponent(self.model_part, VELOCITY_Z_GRADIENT, VORTICITY)

class Pouliot2012EdgeMaterialAccelerationRecoverer(Pouliot2012EdgeGradientRecoverer, recoverer.MaterialAccelerationRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        Pouliot2012EdgeGradientRecoverer.__init__(self, pp, model_part, cplusplus_recovery_tool)
        self.store_full_gradient = self.pp.CFD_DEM.store_full_gradient

    def RecoverMaterialAcceleration(self):
        if self.store_full_gradient:
            self.RecoverGradientOfVelocity()
            self.RecoverMaterialAccelerationFromGradient()
        else:
            self.RecoverGradientOfVelocityComponent(0)
            self.cplusplus_recovery_tool.CalculateVectorMaterialDerivativeComponent(self.model_part, VELOCITY_Z_GRADIENT, ACCELERATION, MATERIAL_ACCELERATION)
            self.RecoverGradientOfVelocityComponent(1)
            self.cplusplus_recovery_tool.CalculateVectorMaterialDerivativeComponent(self.model_part, VELOCITY_Z_GRADIENT, ACCELERATION, MATERIAL_ACCELERATION)
            self.RecoverGradientOfVelocityComponent(2)
            self.cplusplus_recovery_tool.CalculateVectorMaterialDerivativeComponent(self.model_part, VELOCITY_Z_GRADIENT, ACCELERATION, MATERIAL_ACCELERATION)

class Pouliot2012EdgeLaplacianRecoverer(Pouliot2012EdgeMaterialAccelerationRecoverer, recoverer.LaplacianRecoverer):
    def __init__(self, pp, model_part, cplusplus_recovery_tool):
        Pouliot2012EdgeDerivativesRecoverer.__init__(self, pp, model_part, cplusplus_recovery_tool)
        self.element_type = "ComputeVelocityLaplacianSimplex3D"
        self.condition_type = "ComputeLaplacianSimplexCondition3D"
        self.FillUpModelPart(self.element_type)
        self.DOFs = (VELOCITY_LAPLACIAN_X, VELOCITY_LAPLACIAN_Y, VELOCITY_LAPLACIAN_Z)
        self.AddDofs(self.DOFs)

    def RecoverVelocityLaplacian(self):
        print("\nSolving for the laplacian...")
        self.SetToZero(VELOCITY_LAPLACIAN)
        self.recovery_strategy.Solve()
