from abc import abstractmethod
from math import isclose, fabs
import KratosMultiphysics as Kratos

def CreateTerminationCriteria(parameters):
    termination_criteria_type = parameters["termination_criteria"].GetString()

    available_criterians = ["relative_value", "moving_max", "moving_min"]
    if (termination_criteria_type not in available_criterians):
        msg = "Unsupported termination criteria requested. [ \"termination_criteria\" = " + termination_criteria_type
        msg += "\nSupported types: \n\t"
        msg += "\n\t".join(available_criterians)
        raise RuntimeError(msg)

    if (termination_criteria_type == "relative_value"):
        return RelativeToleranceTerminationCriteria(parameters)
    elif (termination_criteria_type == "moving_max"):
        return MovingAverageTerminationCriteria(parameters)
    elif (termination_criteria_type == "moving_min"):
        return MovingMinTerminationCriteria(parameters)

class TerminationCriteria:
    def __init__(self, parameters):
        self.parameters = parameters

    @abstractmethod
    def IsConverged(self, data_logger, optimization_iteration):
        pass

class RelativeToleranceTerminationCriteria(TerminationCriteria):
    def __init__(self, parameters):
        super().__init__(parameters)

        defaults = Kratos.Parameters("""{
            "termination_criteria": "relative_value",
            "relative_tolerance"  : 1e-3
        }""")

        self.parameters.ValidateAndAssignDefaults(defaults)
        self.relative_tolerance = self.parameters["relative_tolerance"].GetDouble()

    def IsConverged(self, data_logger, optimization_iteration):
        relative_change_of_objective_value = data_logger.GetValues("rel_change_objective")[optimization_iteration]
        if abs(relative_change_of_objective_value) < self.relative_tolerance:
            return True

class MovingAverageTerminationCriteria(TerminationCriteria):
    def __init__(self, parameters):
        super().__init__(parameters)

        defaults = Kratos.Parameters("""{
            "termination_criteria"  : "moving_max",
            "moving_max_window"     : 10,
            "list_of_identifiers"   : []
        }""")

        self.parameters.ValidateAndAssignDefaults(defaults)
        self.moving_max_window = self.parameters["moving_max_window"].GetInt()
        self.list_of_identifiers = self.parameters["list_of_identifiers"].GetStringArray()
        self.identifier_max_values = {}

    def IsConverged(self, data_logger, optimization_iteration):
        responses_and_values = data_logger.GetValues("response_value")

        is_converged = optimization_iteration > self.moving_max_window

        if (is_converged):
            Kratos.Logger.Print("")
            for identifier, response_values in responses_and_values.items():
                if (identifier in self.list_of_identifiers):
                    moving_max = response_values[optimization_iteration]
                    for i in range(optimization_iteration - 1, optimization_iteration - self.moving_max_window, -1):
                        moving_max = max(moving_max, response_values[i])

                    if (identifier in self.identifier_max_values.keys()):
                        all_time_max = self.identifier_max_values[identifier]
                        is_converged = is_converged and (all_time_max > moving_max)
                        self.identifier_max_values[identifier] = max(all_time_max, moving_max)
                    else:
                        is_converged = False
                        self.identifier_max_values[identifier] = moving_max

                    Kratos.Logger.PrintInfo("MovingAverageTerminationCriteria", "identifier = {:s}, Current max = {:1.4e}, moving max = {:1.4e}".format(identifier, self.identifier_max_values[identifier], moving_max))

        return is_converged


class MovingMinTerminationCriteria(TerminationCriteria):
    def __init__(self, parameters):
        super().__init__(parameters)

        defaults = Kratos.Parameters("""{
            "termination_criteria"  : "moving_min",
            "moving_min_window"     : 10,
            "list_of_identifiers"   : []
        }""")

        self.parameters.ValidateAndAssignDefaults(defaults)
        self.moving_min_window = self.parameters["moving_min_window"].GetInt()
        self.list_of_identifiers = self.parameters["list_of_identifiers"].GetStringArray()
        self.identifier_min_values = {}

    def IsConverged(self, data_logger, optimization_iteration):
        responses_and_values = data_logger.GetValues("response_value")

        is_converged = optimization_iteration > self.moving_min_window

        if (is_converged):
            Kratos.Logger.Print("")
            for identifier, response_values in responses_and_values.items():
                if (identifier in self.list_of_identifiers):
                    moving_min = response_values[optimization_iteration]
                    for i in range(optimization_iteration - 1, optimization_iteration - self.moving_min_window, -1):
                        moving_min = min(moving_min, response_values[i])

                    if (identifier in self.identifier_min_values.keys()):
                        all_time_min = self.identifier_min_values[identifier]
                        is_converged = is_converged and (all_time_min < moving_min)
                        self.identifier_min_values[identifier] = min(all_time_min, moving_min)
                    else:
                        is_converged = False
                        self.identifier_min_values[identifier] = moving_min

                    Kratos.Logger.PrintInfo("MovingAverageTerminationCriteria", "identifier = {:s}, Current min = {:1.4e}, moving min = {:1.4e}".format(identifier, self.identifier_min_values[identifier], moving_min))

        return is_converged

