from enum import Enum

# Importing the Kratos Library
import KratosMultiphysics as Kratos


def Factory(settings, model):
    if not isinstance(settings, Kratos.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    if not isinstance(model, Kratos.Model):
        raise Exception("Expected input shall be a Model object, encapsulating a json string")

    return ExecuteOperationProcess(model, settings)


class ExecuteOperationProcess(Kratos.Process):
    class __ExecutionPoints(Enum):
        EXECUTE_INITIALIZE = 1
        EXECUTE_BEFORE_SOLUTION_LOOP = 2
        EXECUTE_INITIALIZE_SOLUTION_STEP = 3
        EXECUTE_FINALIZE_SOLUTION_STEP = 4
        EXECUTE_BEFORE_OUTPUT_STEP = 5
        EXECUTE_AFTER_OUTPUT_STEP = 6
        EXECUTE_FINALIZE = 7

    """This process executes given operation at given execution points
    """
    def __init__(self, model: Kratos.Model, params: Kratos.Parameters):
        Kratos.Process.__init__(self)

        params.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.echo_level = params["echo_level"].GetInt()

        operation_name = params["operation_name"].GetString()
        operation_prototype = Kratos.Registry[f"{operation_name}.Prototype"]
        self.operation = operation_prototype.Create(model, params["operation_settings"])

        self.execution_points = []
        allowed_execution_point_names = [data.name.lower() for data in ExecuteOperationProcess.__ExecutionPoints]

        for execution_point_name in params["execution_points"].GetStringArray():
            if execution_point_name not in allowed_execution_point_names:
                raise RuntimeError(f"{execution_point_name} is not valid. Followings are the valid execution points: \n\t" + "\n\t".join(allowed_execution_point_names))

            self.execution_points.append(getattr(ExecuteOperationProcess.__ExecutionPoints, execution_point_name.upper()))

    def GetDefaultParameters(self):
        return Kratos.Parameters("""
                {
                    "help"              : "This process ex",
                    "operation_name"    : "",
                    "operation_settings": {},
                    "execution_points"  : [],
                    "echo_level"        : 0
                }
                """)

    def ExecuteInitialize(self):
        self.__Execute(ExecuteOperationProcess.__ExecutionPoints.EXECUTE_INITIALIZE)

    def ExecuteBeforeSolutionLoop(self):
        self.__Execute(ExecuteOperationProcess.__ExecutionPoints.EXECUTE_BEFORE_SOLUTION_LOOP)

    def ExecuteInitializeSolutionStep(self):
        self.__Execute(ExecuteOperationProcess.__ExecutionPoints.EXECUTE_INITIALIZE_SOLUTION_STEP)

    def ExecuteFinalizeSolutionStep(self):
        self.__Execute(ExecuteOperationProcess.__ExecutionPoints.EXECUTE_FINALIZE_SOLUTION_STEP)

    def ExecuteBeforeOutputStep(self):
        self.__Execute(ExecuteOperationProcess.__ExecutionPoints.EXECUTE_BEFORE_OUTPUT_STEP)

    def ExecuteAfterOutputStep(self):
        self.__Execute(ExecuteOperationProcess.__ExecutionPoints.EXECUTE_AFTER_OUTPUT_STEP)

    def ExecuteFinalize(self):
        self.__Execute(ExecuteOperationProcess.__ExecutionPoints.EXECUTE_FINALIZE)

    def __Execute(self, current_execution_point):
        if current_execution_point in self.execution_points:
            if self.echo_level > 1:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Executing {self.operation.__class__.__name__} at the f{current_execution_point.name.lower()}...")

            self.operation.Execute()

            if self.echo_level > 0:
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Executed {self.operation.__class__.__name__} at the f{current_execution_point.name.lower()}.")

