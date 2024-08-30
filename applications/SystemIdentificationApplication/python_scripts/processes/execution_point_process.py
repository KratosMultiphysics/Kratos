import KratosMultiphysics as Kratos

class ExecutionPointProcess(Kratos.Process):
    def __init__(self, execution_point: str) -> None:
        Kratos.Process.__init__(self)

        self.execution_point = execution_point
        supported_execution_points = ["initialize", "initialize_solution_step", "finalize_solution_step", "finalize"]
        if self.execution_point not in supported_execution_points:
            raise RuntimeError(f"Unsupported execution point = \"{self.execution_point}\". Followings are supported:\n\t" + "\n\t".join(supported_execution_points))

    def Execute(self) -> None:
        raise NotImplementedError("Please implement Execute method in the derived class")

    def ExecuteInitialize(self) -> None:
        if self.execution_point == "initialize":
            self.Execute()

    def ExecuteInitializeSolutionStep(self) -> None:
        if self.execution_point == "initialize_solution_step":
            self.Execute()

    def ExecuteFinalizeSolutionStep(self) -> None:
        if self.execution_point == "finalize_solution_step":
            self.Execute()

    def ExecuteFinalize(self) -> None:
        if self.execution_point == "finalize":
            self.Execute()




