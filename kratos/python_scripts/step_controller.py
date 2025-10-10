from abc import ABC, abstractmethod
import KratosMultiphysics as Kratos

class StepController(ABC):
    """
    Abstract base class for controlling the progression of steps in a simulation.

    Methods
    -------
    Initialize(time_begin: float, time_end: float) -> None
        Prepare the controller for stepping, given the initial and final times.

    GetNextStep(current_time: float, is_converged: bool) -> float
        Determine the next step size based on the current time and convergence status of the current step.

    IsCompleted(current_time: float, is_converged: bool) -> bool
        Check if the stepping process has been completed.
    """
    @abstractmethod
    def Initialize(self, time_begin: float, time_end: float) -> None:
        """
        Initializes the step controller with the specified start and end times for current time step from the solver.

        Args:
            time_begin (float): The starting time of the current time step from the solver.
            time_end (float): The ending time of the current time step from the solver.

        Returns:
            None
        """
        pass

    @abstractmethod
    def GetNextStep(self, current_time: float, is_converged: bool) -> float:
        """
        Determines and returns the next step based on the current simulation time and convergence status.

        Args:
            current_time (float): The current time in the simulation.
            is_converged (bool): Flag indicating whether the previous step has converged.

        Returns:
            float: The time for the next step.
        """
        pass

    @abstractmethod
    def IsCompleted(self, current_time: float, is_converged: bool) -> bool:
        """
        Called to check if the stepping process is completed.

        Args:
            current_time (float): The current simulation time.
            is_converged (bool): Indicates whether the solution has converged at the current step.

        Returns:
            bool: True if the stepping is completed, otherwise false.
        """
        pass

class DefaultStepController(StepController):
    """
    DefaultStepController is a step controller implementation for Kratos that always returns the final time as the next step,
    effectively disabling intermediate stepping.

    Args:
        parameters (Kratos.Parameters): Configuration parameters for the step controller.

    Methods:
        Initialize(start_time: float, time_end: float) -> None:
            Initializes the controller with the end time of the simulation.

        GetNextStep(current_time: float, is_converged: bool) -> float:
            Returns the end time as the next step, regardless of the current time or convergence status.

        IsCompleted(current_time: float, is_converged: bool) -> bool:
            Always returns True, indicating that stepping is completed.
    """
    def __init__(self, parameters: Kratos.Parameters) -> None:
        default_parameters = Kratos.Parameters("""{
            "type": "default_step_controller"
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

    def Initialize(self, _: float, time_end: float) -> None:
        self.__time_end = time_end

    def GetNextStep(self, _: float, __: bool) -> float:
        return self.__time_end

    def IsCompleted(self, _: float, __: bool) -> bool:
        return True

class GeometricStepController(StepController):
    """
    GeometricStepController implements an adaptive time-stepping controller based on geometric progression for stepping in simulations.

    This controller adjusts the time increment (`delta_time`) dynamically based on the convergence or divergence of the solution at each step. It increases the time step after a specified number of successful convergences and decreases it after a failed attempt, within user-defined bounds.

    Class Methods:
        GetDefaultParameters() -> Kratos.Parameters:
            Returns the default parameters for the controller, including factors for convergence/divergence, initial/min/max time steps, and thresholds for incrementing or terminating the step size.

    Constructor:
        __init__(settings: Kratos.Parameters):
            Initializes the controller with user-provided or default settings, validating and storing parameters for each interval.

    Methods:
        Initialize(time_begin: float, time_end: float) -> None:
            Prepares the controller for a new stepping process, resetting counters and setting initial parameters.

        GetNextStep(current_time: float, is_converged: bool) -> float:
            Determines the next time step based on the convergence status of the current time step. Increments or decrements the time step according to the configured factors and thresholds.

        IsCompleted(current_time: float, is_converged: bool) -> bool:
            Checks if the stepping process is completed, i.e., the solution has converged and the end time is reached.

        __SetSubSteppingParameters() -> None:
            Selects and sets the appropriate sub-stepping parameters based on the current interval.

    Attributes:
        __divergence_factor: Factor by which the time step is reduced after a failed attempt.
        __convergence_factor: Factor by which the time step is increased after successful attempts.
        __delta_t_init: Initial time step size.
        __delta_t_min: Minimum allowed time step size.
        __delta_t_max: Maximum allowed time step size.
        __max_number_of_sub_steps: Maximum number of sub-steps allowed.
        __number_of_successful_attempts_for_increment: Number of successful steps required to increase the time step.
        __number_of_failed_attempts_for_termination: Number of failed steps allowed before termination.
        __sub_stepping_iteration: Counter for the current sub-step iteration.
        __success_full_attempts_count: Counter for successful attempts.
        __failed_attempts_count: Counter for failed attempts.
        __delta_time: Current time step size.
        __time_begin: Start time of the stepping process.
        __time_end: End time of the stepping process.

    Raises:
        RuntimeError: If the maximum number of sub-steps or failed attempts is exceeded, or if the minimum allowed time step is reached.
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

    def __init__(self, settings: Kratos.Parameters) -> None:
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

    def Initialize(self, time_begin: float, time_end: float) -> None:
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