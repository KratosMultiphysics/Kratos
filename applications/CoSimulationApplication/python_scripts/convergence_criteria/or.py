from KratosMultiphysics.CoSimulationApplication.convergence_criteria.combined import ConvergenceCriterionCombined
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return ConvergenceCriterionOr(parameters)


class ConvergenceCriterionOr(ConvergenceCriterionCombined):
    def IsSatisfied(self):
        is_satisfied = False
        for convergence_criterion in self.convergence_criteria:
            is_satisfied = is_satisfied or convergence_criterion.IsSatisfied()

        return is_satisfied
