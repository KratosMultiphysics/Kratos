from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos

class StepController(ABC):
    @abstractmethod
    def Initialize(self, time_begin: float, time_end: float) -> None:
        pass

    @abstractmethod
    def GetSubStep(self, current_time: float, is_converged: bool) -> float:
        pass

    @abstractmethod
    def SubSteppingCompleted(self, current_time: float, is_converged: bool) -> bool:
        pass

class DefaultStepController(StepController):
    def __init__(self, parameters: Kratos.Parameters) -> None:
        default_parameters = Kratos.Parameters("""{
            "type": "default_step_controller"
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

    def Initialize(self, _: float, time_end: float) -> None:
        self.__time_end = time_end

    def GetSubStep(self, _: float, __: bool) -> float:
        return self.__time_end

    def SubSteppingCompleted(self, _: float, __: bool) -> bool:
        return True

class GeometricStepController(StepController):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"               : "geometric_step_controller",
            "list_of_controllers": [
                {
                    "interval"                                   : [0.0, "End"],
                    "divergence_factor"                          : 0.25,
                    "convergence_factor"                         : 1.5,
                    "delta_t_init"                               : 1.0,
                    "delta_t_min"                                : 1.0,
                    "delta_t_max"                                : 1.0,
                    "max_number_of_sub_steps"                    : 100,
                    "number_of_successful_attempts_for_increment": 2,
                    "number_of_failed_attempts_for_termination"  : 5

                }
            ]
        }""")

    def __init__(self, settings: Kratos.Parameters) -> None:
        default_parameters = self.GetDefaultParameters()

        settings.ValidateAndAssignDefaults(default_parameters)

        self.__list_of_step_controllers: 'list[tuple[Kratos.IntervalUtility, Kratos.Parameters]]' = []
        for controller_settings in settings["list_of_controllers"].values():
            controller_settings.ValidateAndAssignDefaults(default_parameters["list_of_controllers"].values()[0])
            self.__list_of_step_controllers.append((Kratos.IntervalUtility(controller_settings), controller_settings.Clone()))

    def Initialize(self, time_begin: float, time_end: float) -> None:
        self.__time_begin = time_begin
        self.__time_end = time_end

        self.__SetSubSteppingParameters()

        self.__sub_stepping_iteration = 0
        self.__success_full_attempts_count = 0
        self.__failed_attempts_count = 0
        self.__delta_time = self.__delta_t_init

    def GetSubStep(self, current_time: float, is_converged: bool) -> float:

        self.__sub_stepping_iteration += 1
        if self.__sub_stepping_iteration >= self.__max_number_of_sub_steps:
            raise RuntimeError(f"{self.__class__.__name__}: Reached maximum number of sub steps [ max number of sub steps = {self.__max_number_of_sub_steps} ].")

        if is_converged:
            # the sub-step is converged. So try the next sub step
            self.__failed_attempts_count = 0
            self.__success_full_attempts_count += 1

            if self.__success_full_attempts_count >= self.__number_of_successful_attempts_for_increment:
                # now we increment the step
                self.__delta_time = min(self.__convergence_factor * self.__delta_time, self.__delta_t_max)
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Incrementing the delta time to {self.__delta_time}")

            return min(self.__time_end, current_time + self.__delta_time)
        else:
            # the sub-step is not converged.
            self.__failed_attempts_count += 1
            self.__success_full_attempts_count = 0

            if (self.__failed_attempts_count >= self.__number_of_failed_attempts_for_termination):
                raise RuntimeError(f"Reached maximum number of failed attempts [ {self.__failed_attempts_count} / {self.__number_of_failed_attempts_for_termination} ]")
            else:
                self.__delta_time = self.__divergence_factor * self.__delta_time
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Decrementing the delta time to {self.__delta_time}")
                if self.__delta_time <= self.__delta_t_min:
                    raise RuntimeError(f"Reached the minimum allowed delta time [ delta_time = {self.__delta_time}, allowed minimum delta_time = {self.__delta_t_min} ]")

            new_time = current_time + self.__delta_time
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Solving sub step {self.__sub_stepping_iteration} for time = {new_time:0.6e}")
            return new_time

    def SubSteppingCompleted(self, current_time: float, is_converged: bool) -> bool:
        return is_converged and abs(current_time - self.__time_end) <= (self.__time_end - self.__time_begin) * 1e-9

    def __SetSubSteppingParameters(self) -> None:
        for interval_utility, controller_settings in self.__list_of_step_controllers:
            if interval_utility.IsInInterval(self.__time_end):
                    self.__divergence_factor = controller_settings["divergence_factor"].GetDouble()
                    self.__convergence_factor = controller_settings["convergence_factor"].GetDouble()
                    self.__delta_t_init = controller_settings["delta_t_init"].GetDouble()
                    self.__delta_t_min = controller_settings["delta_t_min"].GetDouble()
                    self.__delta_t_max = controller_settings["delta_t_max"].GetDouble()
                    self.__max_number_of_sub_steps = controller_settings["max_number_of_sub_steps"].GetInt()
                    self.__number_of_successful_attempts_for_increment = controller_settings["number_of_successful_attempts_for_increment"].GetInt()
                    self.__number_of_failed_attempts_for_termination = controller_settings["number_of_failed_attempts_for_termination"].GetInt()
                    return None

        raise RuntimeError(f"")

def Factory(parameters: Kratos.Parameters) -> StepController:
    if not parameters.Has("type"):
        parameters.AddString("type", "default_step_controller")

    step_controller_map = {
        "default_step_controller": DefaultStepController,
        "geometric_step_controller": GeometricStepController
    }

    controller_type = parameters["type"].GetString()
    if controller_type not in step_controller_map.keys():
        raise RuntimeError(f"Unsupported controller type = \"{controller_type}\". Followings are supported:\n\t" + "\n\t".join(step_controller_map.keys()))

    return step_controller_map[controller_type](parameters)