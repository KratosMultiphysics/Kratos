from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos

class StepController(ABC):
    """\
    @brief Abstract base class for controlling the progression of steps in a simulation.

    @details
    Defines the interface used by all step controllers to manage time stepping
    behaviour during a simulation. Implementations must provide logic for
    determining the next step size, and checking completion.
    """

    @abstractmethod
    def GetNextStep(self, current_time: float, is_converged: bool) -> float:
        """\
        @brief Compute the next stepping time.

        @param[in] current_time The current simulation time.
        @param[in] is_converged True if the previous step converged successfully.

        @return The proposed time for the next step.
        """
        pass

    @abstractmethod
    def IsCompleted(self, current_time: float, is_converged: bool) -> bool:
        """\
        @brief Query whether the stepping procedure has finished.

        @param[in] current_time The current simulation time.
        @param[in] is_converged Flag indicating convergence at the current step.

        @return True when no further steps are required, false otherwise.
        """
        pass

class DefaultStepController(StepController):
    """\
    @brief A trivial controller that jumps directly to the final time.

    @details
    This step controller ignores intermediate steps and always returns the
    final simulation time as the next step. It is useful when no sub-stepping
    behaviour is desired.
    """

    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type": "default_step_controller"
        }""")

    def __init__(self, _: float, time_end: float, parameters: Kratos.Parameters) -> None:
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.__time_end = time_end

    def GetNextStep(self, _: float, __: bool) -> float:
        return self.__time_end

    def IsCompleted(self, _: float, __: bool) -> bool:
        return True

class GeometricStepController(StepController):
    """\
    @brief Adaptive geometric time-step controller.

    @details
    A controller that modifies the time increment according to a geometric
    progression. When consecutive steps converge, the step size is multiplied
    by a convergence factor up to a maximum. Conversely, on divergence the step
    size is scaled down by a divergence factor until a minimum is reached. The
    controller also tracks the number of sub-steps and failed attempts, throwing
    runtime errors if limits are exceeded.
    """
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"                                       : "geometric_step_controller",
            "divergence_factor"                          : 0.25,
            "convergence_factor"                         : 1.5,
            "delta_t_init"                               : 1.0,
            "delta_t_min"                                : 1.0,
            "delta_t_max"                                : 1.0,
            "max_number_of_sub_steps"                    : 100,
            "number_of_successful_attempts_for_increment": 2,
            "number_of_failed_attempts_for_termination"  : 5
        }""")

    def __init__(self, time_begin: float, time_end: float, settings: Kratos.Parameters) -> None:
        default_parameters = self.GetDefaultParameters()

        settings.ValidateAndAssignDefaults(default_parameters)

        self.__divergence_factor = settings["divergence_factor"].GetDouble()
        self.__convergence_factor = settings["convergence_factor"].GetDouble()
        self.__delta_t_init = settings["delta_t_init"].GetDouble()
        self.__delta_t_min = settings["delta_t_min"].GetDouble()
        self.__delta_t_max = settings["delta_t_max"].GetDouble()
        self.__max_number_of_sub_steps = settings["max_number_of_sub_steps"].GetInt()
        self.__number_of_successful_attempts_for_increment = settings["number_of_successful_attempts_for_increment"].GetInt()
        self.__number_of_failed_attempts_for_termination = settings["number_of_failed_attempts_for_termination"].GetInt()

        self.__time_begin = time_begin
        self.__time_end = time_end

        self.__sub_stepping_iteration = 0
        self.__success_full_attempts_count = 0
        self.__failed_attempts_count = 0
        self.__delta_time = self.__delta_t_init

    def GetNextStep(self, current_time: float, is_converged: bool) -> float:

        self.__sub_stepping_iteration += 1
        if self.__sub_stepping_iteration >= self.__max_number_of_sub_steps:
            raise RuntimeError(f"{self.__class__.__name__}: Reached maximum number of Load steps [ max number of Load steps = {self.__max_number_of_sub_steps} ].")

        if is_converged:
            # the sub-step is converged. So try the next Step
            self.__failed_attempts_count = 0
            self.__success_full_attempts_count += 1

            if self.__success_full_attempts_count > self.__number_of_successful_attempts_for_increment:
                # now we increment the step
                self.__delta_time = min(self.__convergence_factor * self.__delta_time, self.__delta_t_max)
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Incrementing the delta time to {self.__delta_time}")

            return min(self.__time_end, current_time + self.__delta_time)
        else:
            # the sub-step is not converged.
            self.__failed_attempts_count += 1
            self.__success_full_attempts_count = 0

            if (self.__failed_attempts_count > self.__number_of_failed_attempts_for_termination):
                raise RuntimeError(f"Reached maximum number of failed attempts [ {self.__failed_attempts_count} / {self.__number_of_failed_attempts_for_termination} ]")
            else:
                self.__delta_time = self.__divergence_factor * self.__delta_time
                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Decrementing the delta time to {self.__delta_time}")
                if self.__delta_time <= self.__delta_t_min:
                    raise RuntimeError(f"Reached the minimum allowed delta time [ delta_time = {self.__delta_time}, allowed minimum delta_time = {self.__delta_t_min} ]")

            new_time = current_time + self.__delta_time
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Solving Step {self.__sub_stepping_iteration} for time = {new_time:0.6e}")
            return new_time

    def IsCompleted(self, current_time: float, is_converged: bool) -> bool:
        return is_converged and abs(current_time - self.__time_end) <= (self.__time_end - self.__time_begin) * 1e-9

class CompoundStepController(StepController):
    """\
    @brief CompoundStepController handles different step controllers over different time intervals.

    @details
    This class allows for adaptive time stepping strategies that change based on the current
    time period. It selects the appropriate step controller from a list based on which interval
    the current time falls within.

    @throws RuntimeError If no step controller is found for the given time period.
    """
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "type"                    : "compound_step_controller",
            "list_of_step_controllers": [
                {
                    "interval": [0.0, "End"],
                    "settings": {}
                }
            ]
        }""")

    def __init__(self, time_begin: float, time_end: float, settings: Kratos.Parameters) -> None:
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.step_controller = None
        for step_controller_setting in reversed(settings["list_of_step_controllers"].values()):
            current_settings: Kratos.Parameters = step_controller_setting.Clone()
            current_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["list_of_step_controllers"].values()[0])
            if Kratos.IntervalUtility(current_settings).IsInInterval(time_end):
                self.step_controller = Factory(time_begin, time_end, current_settings["settings"])
                break

        if self.step_controller is None:
            raise RuntimeError(f"No step controller is found for the time period [{time_begin}, {time_end}].")

    def GetNextStep(self, current_time: float, is_converged: bool) -> float:
        return self.step_controller.GetNextStep(current_time, is_converged)

    def IsCompleted(self, current_time: float, is_converged: bool) -> bool:
        return self.step_controller.IsCompleted(current_time, is_converged)

def Factory(time_begin: float, time_end: float, parameters: Kratos.Parameters) -> StepController:
    if not parameters.Has("type"):
        parameters.AddString("type", "default_step_controller")

    step_controller_map = {
        "default_step_controller": DefaultStepController,
        "geometric_step_controller": GeometricStepController,
        "compound_step_controller": CompoundStepController
    }

    controller_type = parameters["type"].GetString()
    if controller_type not in step_controller_map.keys():
        raise RuntimeError(f"Unsupported controller type = \"{controller_type}\". Followings are supported:\n\t" + "\n\t".join(step_controller_map.keys()))

    return step_controller_map[controller_type](time_begin, time_end, parameters)