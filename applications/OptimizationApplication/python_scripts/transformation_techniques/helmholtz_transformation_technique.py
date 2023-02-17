import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.transformation_techniques.transformation_technique import TransformationTechnique
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import ContainerVariableDataHolderUnion
from KratosMultiphysics.OptimizationApplication.execution_policies.stepping_analysis_execution_policy import SteppingAnalysisExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.helmholtz_analysis import HelmholtzAnalysis
from KratosMultiphysics.OptimizationApplication.solvers.helmholtz_solver import HelmholtzSolver

class HelmholtzTransformationTechnique(TransformationTechnique):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__()

        default_parameters = Kratos.Parameters("""{
            "helmholtz_execution_policy_name": "",
            "filter_radius"                  : 1.0,
            "fixed_model_part_names"         : []
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)
        self.model = model
        self.filter_radius = parameters["filter_radius"].GetDouble()
        self.helmholtz_execution_policy_decorator: ExecutionPolicyDecorator = optimization_info.GetOptimizationProcess(ExecutionPolicyDecorator, parameters["helmholtz_execution_policy_name"].GetString())
        self.fixed_model_part_names = parameters["fixed_model_part_names"].GetStringArray()

        execution_policy: SteppingAnalysisExecutionPolicy = self.helmholtz_execution_policy_decorator.GetExecutionPolicy()
        if not isinstance(execution_policy, SteppingAnalysisExecutionPolicy):
            raise RuntimeError(f"Helmholtz tranformation requires the Helmholtz analysis to be in a SteppingAnalysisExecutionPolicy. The provided execution policy \"{self.helmholtz_execution_policy_decorator.GetExecutionPolicyName()}\" is of type {execution_policy.__class__.__name__}.")

        self.helmholtz_analysis: HelmholtzAnalysis = execution_policy.analysis
        if not isinstance(self.helmholtz_analysis, HelmholtzAnalysis):
            raise RuntimeError(f"Helmholtz tranformation requires the SteppingAnalysisExecutionPolicy to have a Helmholtz analysis. The provided analysis is of type \"{self.helmholtz_execution_policy_decorator.GetExecutionPolicyName()}\" is of type {self.helmholtz_analysis.__class__.__name__}.")

        self.helmholtz_solver: HelmholtzSolver = self.helmholtz_analysis._GetSolver()

    def TransformSensitivity(self, container_variable_data_holder: ContainerVariableDataHolderUnion):
        # set filter radius
        radius_container = KratosOA.ElementPropertiesContainerVariableDataHolder(self.helmholtz_solver.GetComputingModelPart())
        radius_container.SetDataForContainerVariableToZero(KratosOA.HELMHOLTZ_RADIUS_DENSITY)
        radius_container += self.filter_radius
        radius_container.AssignDataToContainerVariable(KratosOA.HELMHOLTZ_RADIUS_DENSITY)

        self.helmholtz_solver.SetData(container_variable_data_holder)
        self.helmholtz_execution_policy_decorator.Execute()
        self.helmholtz_solver.GetData(container_variable_data_holder)

    def TransformUpdate(self, container_variable_data_holder: ContainerVariableDataHolderUnion):
        # set filter radius
        radius_container = KratosOA.ElementPropertiesContainerVariableDataHolder(self.helmholtz_solver.GetComputingModelPart())
        radius_container.SetDataForContainerVariableToZero(KratosOA.HELMHOLTZ_RADIUS_DENSITY)
        radius_container += self.filter_radius
        radius_container.AssignDataToContainerVariable(KratosOA.HELMHOLTZ_RADIUS_DENSITY)

        container_model_part: Kratos.ModelPart = container_variable_data_holder.GetModelPart()

        nodal_values = KratosOA.HistoricalContainerVariableDataHolder(container_model_part)
        self.helmholtz_solver.ConvertSolverContainer(container_variable_data_holder, nodal_values)

        mass_conribution = KratosOA.HistoricalContainerVariableDataHolder(container_model_part)
        KratosOA.ContainerVariableDataHolderUtils.ComputeVariableDataHolderProductWithEntityMatrix(mass_conribution, nodal_values, KratosOA.HELMHOLTZ_MASS_MATRIX, self.helmholtz_solver.GetComputingModelPart().Elements)

        self.helmholtz_solver.SetData(mass_conribution)

        for fixed_model_part_name in self.fixed_model_part_names:
            fixed_model_part = self.model[fixed_model_part_name]
            fixed_values = KratosOA.HistoricalContainerVariableDataHolder(fixed_model_part)
            fixed_values.SetDataForContainerVariableToZero(self.helmholtz_solver.GetDofVariable())
            self.helmholtz_solver.FixDofs(fixed_values)

        self.helmholtz_execution_policy_decorator.Execute()

        self.helmholtz_solver.GetData(container_variable_data_holder)

