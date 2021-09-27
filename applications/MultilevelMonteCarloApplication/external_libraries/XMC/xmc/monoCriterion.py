from warnings import warn
from typing import Callable, Optional, List

# XMC imports
from xmc.tools import dynamicImport
from xmc.methodDefs_monoCriterion import criterionFunctions


class MonoCriterion:
    """
    A generic, elementary criterion object.
    """

    _criterion: Callable
    """Function to perform the comparison. Takes two or three floats as parameters. Returns boolean."""

    _tolerance: Optional[float]
    """Tolerance to be used for comparison. If ``None``, the :py:meth:`~.flag` method expected two input arguments."""

    def __init__(self, criterionFunction: str, tolerance: Optional[float] = None):
        """
        - criterionFunction is a necessary input argument. It is a string of the full name of the criterion function, e.g. xmc.methodDefs_monoCriterion.criterionFunctions.isLowerThanOrEqualTo.
        - tolerance is not always necessary; it depends on the chosen criterionFunction.
        """
        self._criterion = dynamicImport(criterionFunction)
        if tolerance is None:
            return

        if isinstance(tolerance, float):
            self._tolerance = tolerance
        # What follow is to handle incorrect formats
        elif isinstance(tolerance, list) and all(isinstance(e, float) for e in tolerance):
            if len(tolerance) == 1:
                # This is the former format: [tolerance_as_single_float]
                # Accept it, but warn that it is deprecated
                warn(
                    "Passing single float in list is deprecated for 'tolerance'. Pass float directly.",
                    FutureWarning,
                )
                self._tolerance = tolerance[0]
            else:
                # List of multiple tolerances. Not supported by this class.
                raise NotImplementedError(
                    (
                        "MonoCriterion does not support a sequence of tolerances. "
                        "Consider using CriterionSequence."
                    )
                )
        else:
            # Unexpected. Entirely wrong type.
            raise TypeError(
                f"Expected parameter tolerance to be float; received {type(tolerance)}"
            )

    # TODO - Monocriterion.flag should be able to run a criterion
    # without a tolerance. Hence, b should not default to None and
    # the cases should be appropriately changed
    def flag(self, a: float, b: float = None) -> bool:
        """
        Takes at most two input arguments (a and b) and re-
        turns a boolean. Calls criterion(a, b, tolerance).
        Some criterion functions may need only two input arguments, e.g. (a,tolerance) -> a<tolerance.
        """
        if b is None:
            return self._criterion(a, self.tolerance)
        else:
            return self._criterion(a, b, self.tolerance)

    @property
    def tolerance(self) -> float:
        """The tolerance used by the criterion"""
        return self._tolerance


class CriterionSequence(MonoCriterion):
    """MonoCriterion with a sequence of tolerances; can be updated."""

    _tolerances: List[float]
    """Sequence of tolerances, in order in which to consider them."""

    def __init__(self, criterionFunction: str, tolerances: List[float]) -> None:
        """
        Create CriterionSequence object.

        :param criterionFunction: import path for attribute :py:attr:`.~_criterion`.
        :param tolerances: list of tolerances, for attribute :py:attr:`.~_tolerances`.
        """

        super().__init__(criterionFunction)
        self._tolerances = tolerances

    @property
    def tolerance(self) -> float:
        """The tolerance used by the criterion"""
        return self._tolerances[0]

    def updateTolerance(self) -> None:
        """Starts using the next tolerance in sequence. Does nothing if there is no next tolerance"""
        if len(self._tolerances) > 1:
            del self._tolerances[0]
