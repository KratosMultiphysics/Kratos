import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_constraint import StandardizedConstraint
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.LinearSolversApplication.dense_linear_solver_factory import ConstructSolver


class GradientProjection(PythonSolver):
    """
        A classical gradient projection algorithm to solve unconstrainted optimization problems.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "model_part_name"   : "OptimizationModelPart",
            "objectives"        : [],
            "constraints"       : [],
            "controls"          : [],
            "echo_level"        : 0,
            "settings"          : {
                "gradient_scaling"       : "inf_norm",
                "echo_level"             : 0,
                "step_size"              : float, 
                "max_iter"               : float,
                "linear_solver_settings" : {},
            }
        }""")
    
    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(model, parameters)

        self.model = model
        self.parameters = parameters
        self.optimization_problem = optimization_problem

        self.master_control = MasterControl() # Need to fill it with controls
        control_list = parameters["controls"]

        for control in control_list:
            control = Control() # Use Control Factory 
            self.master_control.AddControl(control)

        algorithm_parameters = parameters["settings"]
        algorithm_parameters.ValidateAndAssignDefaults(self.GetDefaultParameters()["settings"])

        self.gradient_scaling = algorithm_parameters["gradient_scaling"].GetString()
        self.echo_level = algorithm_parameters["echo_level"].GetInt()

        self.step_size = algorithm_parameters["step_size"].GetInt()
        self.__max_iter = algorithm_parameters["max_iter"].GetInt()

        self.__objective = StandardizedObjective(parameters["objectives"][0], self.master_control, optimization_problem)
        self.__constraints_list = []
        for constraint_param in parameters["constraints"]:
            constraint = StandardizedConstraint(constraint_param, self.master_control, optimization_problem)
            self.__constraints_list.append(constraint)
        self.control_field = None

        default_linear_solver_settings = Kratos.Parameters("""{
            "solver_type": "LinearSolversApplication.dense_col_piv_householder_qr"
        }""")
        algorithm_parameters["linear_solver_settings"].ValidateAndAssignDefaults(default_linear_solver_settings)
        self.linear_solver = ConstructSolver(algorithm_parameters["linear_solver_settings"])

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self):
        if len(self.GetObjectives()) > 1:
            raise RuntimeError(f"{self.__class__.__name__} only supports single objective optimizations.")
        
    def Initialize(self):
        super().Initialize()
        self.master_control.Initialize()

        self.converged = False
        self.control_field = None
        
        self.control_field = self.master_control.GetEmptyControlFields() # GetInitialControlFields() later

    def Finalize(self):
        self.opt_iter = 1
        return super().Finalize()
    
    def ComputeSearchDirection(self, obj_grad, N, g_a) -> KratosOA.ContainerExpression.CollectiveExpressions:
        number_of_active_constraints = len(N)

        # If no active constraints, then do steepest descent
        if number_of_active_constraints == 0:
            return obj_grad * -1.0 

        # Comput the projected search direction and correction
        # 1. N transpose N
        NTN = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints)
        for i in range(number_of_active_constraints):
            for j in range(i, number_of_active_constraints):
                NTN[i, j] = KratosOA.ContainerExpressionUtils.InnerProduct(N[i] * N[j])
                NTN[j, i] = NTN[i, j]
        # 2. inverse NTN
        NTN_inverse = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints)
        I = Kratos.Matrix(number_of_active_constraints, number_of_active_constraints)
        for i in range(number_of_active_constraints):
            I[i, i] = 1.0
        self.linear_solver.Solve(NTN, NTN_inverse, I)

        NT_grad_obj = Kratos.Vector(number_of_active_constraints)
        for i in range(number_of_active_constraints):
            NT_grad_obj[i] = KratosOA.ContainerExpressionUtils.InnerProduct(N[i], obj_grad)

        Lambda = Kratos.Vector(number_of_active_constraints)
        Lambda = NTN_inverse * NT_grad_obj
        
        NLambda = KratosOA.ContainerExpression.CollectiveExpressions()
        NLambda.SetToZero()
        for i in range(number_of_active_constraints):
            NLambda += Lambda[i] * N[i]


        projection_direction = - (obj_grad + self.__CollectiveListVectorProduct(N, NTN_inverse * self.__CollectiveListCollectiveProduct(N, obj_grad)))
        correction_move = - self.__CollectiveListVectorProduct(N, NTN_inverse * g_a)

        return projection_direction, correction_move
    
    def LineSearch(self, search_direction) -> float:
        return self.step_size / KratosOA.ContainerExpressionUtils.NormInf(search_direction)

    def SolveSolutionStep(self) -> bool:
        self.Initialize()

        while not self.converged:

            obj_val = self.__objective.CalculateStandardizedValue(self.control_field) # obj_val is typically not used. It is needed if we use Line Search Technique
            obj_grad = self.__objective.CalculateStandardizedGradient()

            N: 'list[KratosOA.ContainerExpression.CollectiveExpressions]' = []
            g_a = []
            for constraint in self.__constraints_list:
                constr_value = constraint.CalculateStandardizedValue(self.control_field)
                if constr_value >= 0.0:
                    constr_grad = self.__objective.CalculateStandardizedGradient()
                    if KratosOA.ContainerExpressionUtils.NormInf(constr_grad) > 1e-6:
                        N.append(constr_grad)
                        g_a.append(constr_value)

            projection_direction, correction_move = self.ComputeSearchDirection(obj_grad, N, g_a)
            alpha = self.LineSearch(projection_direction)
            self.control_field += projection_direction * alpha + correction_move # we can add limitter for correction move to linit the update or so

            self.converged = self.CheckConvergence()

        self.Finalize()

    def CheckConvergence(self) -> bool:
        return True if self.opt_iter >= self.__max_iter else False

    def __CollectiveListCollectiveProduct(self, collective_list: 'list[KratosOA.CollectiveExpressions]', other_collective: KratosOA.CollectiveVariableDataHolder) -> Kratos.Vector:
        result = Kratos.Vector(len(collective_list))

        for i, collective_list_item in enumerate(collective_list):
            result[i] = KratosOA.ContainerVariableDataHolderUtils.InnerProduct(collective_list_item, other_collective)
        return result

    def __CollectiveListVectorProduct(self, collective_list: 'list[KratosOA.CollectiveExpressions]', vector: Kratos.Vector) -> KratosOA.CollectiveExpressions:
        if len(collective_list) != vector.Size():
            raise RuntimeError(f"Collective list size and vector size mismatch. [ Collective list size = {len(collective_list)}, vector size = {vector.Size()}]")
        if len(collective_list) == 0:
            raise RuntimeError(f"Collective lists cannot be empty.")

        result = collective_list[0].CloneWithDataInitializedToZero()
        for i, collective_list_item in enumerate(collective_list):
            result += collective_list_item * vector[i]

        return result