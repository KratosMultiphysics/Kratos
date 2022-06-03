# XMC imports
from xmc.tools import dynamicImport
from xmc.methodDefs_monoCriterion import criterionFunctions


class MonoCriterion:
    """
    A generic, elementary criterion object.
    """

    def __init__(self, criterion_function, tolerance_values=None):
        """
        - criterion_function is a necessary input argument. It is a string of the full name of the criterion function, e.g. xmc.methodDefs_monoCriterion.criterionFunctions.isLowerThanOrEqualTo.
        - tolerance_value is not always necessary; it depends on the chosen criterion_function.
        """
        self._criterion = dynamicImport(criterion_function)

        if len(tolerance_values) == 0:
            self.tolerance = tolerance_values
        elif len(tolerance_values) == 1:
            self.tolerance = tolerance_values[0]
        else:
            self.tolerance = tolerance_values.pop(0)
        self.tolerances = tolerance_values

    # TODO - Monocriterion.flag should be able to run a criterion
    # without a tolerance. Hence, b should not default to None and
    # the cases should be appropriately changed
    def flag(self, a, b=None):
        """
        Takes at most two input arguments (a and b) and re-
        turns a boolean. Calls criterion(a, b, tolerance).
        Some criterion functions may need only two input arguments, e.g. (a,tolerance) -> a<tolerance.
        """
        if b is None:
            return self._criterion(a, self.tolerance)
        else:
            return self._criterion(a, b, self.tolerance)
