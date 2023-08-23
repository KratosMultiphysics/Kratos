#    |  /           |
#    ' /   __| _` | __|  _ \   __|
#    . \  |   (   | |   (   |\__ `
#   _|\_\_|  \__,_|\__|\___/ ____/
#                   Multi-Physics
#
#  License:		 BSD License
#					 license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#                   Suneth Warnakulasuriya
#

# additional imports
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger

# ==============================================================================
def CreateController(execution_policies_settings,model,model_parts_controller):
    return AnalysesController(execution_policies_settings,model,model_parts_controller)

# ==============================================================================
class AnalysesController:
    # --------------------------------------------------------------------------
    def __init__(self, execution_policies_settings, model, model_parts_controller):
        self.execution_policies_settings = execution_policies_settings
        self.model_parts_controller = model_parts_controller
        self.model = model

        self.execution_policy_decorators = {}
        for execution_policy_settings in self.execution_policies_settings:
            execution_policy_wrapper = ExecutionPolicyDecorator(self.model, execution_policy_settings, OptimizationProblem())
            if not execution_policy_wrapper.GetName() in self.execution_policy_decorators.keys():
                self.execution_policy_decorators[execution_policy_wrapper.GetName()] = execution_policy_wrapper
            else:
                raise RuntimeError(f"Found already existing execution policy with the name \"{execution_policy_wrapper.GetExecutionPolicyName()}\". Please provide unique names.")

    # --------------------------------------------------------------------------
    def GetExecutionPolicyDecorator(self, execution_policy_name: str) -> ExecutionPolicyDecorator:
        if not execution_policy_name in self.execution_policy_decorators.keys():
            raise RuntimeError("AnalysesController: Try to get an execution policy {} which does not exist.".format(execution_policy_name))
        else:
            return self.execution_policy_decorators[execution_policy_name]

    # --------------------------------------------------------------------------
    def GetAnalysis(self, execution_policy_name: str):
        return self.GetExecutionPolicyDecorator(execution_policy_name).GetExecutionPolicy().analysis

    # --------------------------------------------------------------------------
    def Initialize(self):
        for execution_policy_wrapper in self.execution_policy_decorators.values():
            execution_policy_name = execution_policy_wrapper.GetName()
            with TimeLogger(self.__class__.__name__, f"Initializing {execution_policy_name}...", f"Finished initializing {execution_policy_name}"):
                execution_policy_wrapper.Initialize()

    # --------------------------------------------------------------------------
    def RunAnalysis(self, execution_policy_name: str):
        with TimeLogger(self.__class__.__name__, f"Starting {execution_policy_name}...", f"Finished execution of {execution_policy_name}."):
            self.GetExecutionPolicyDecorator(execution_policy_name).Execute()

    # --------------------------------------------------------------------------
    def RunAll(self):
        for name in self.analyses.keys():
            self.RunAnalysis(name)

    # --------------------------------------------------------------------------
    def RunAnalyses(self, analyses_name):
        if not isinstance(analyses_name, list):
            raise RuntimeError("AnalysesController: RunAnalyses requires list of analysis names")
        for analysis_name in analyses_name:
            self.RunAnalysis(analysis_name)



