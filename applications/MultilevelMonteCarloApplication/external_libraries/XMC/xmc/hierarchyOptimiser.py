from typing import Union, Callable, List

from exaqute import *

# Import XMC methods
from xmc.tools import dynamicImport, mergeTwoListsIntoOne
from xmc.methodDefs_hierarchyOptimiser import (
    optimalSampleNumbers,
    optimalIndexSet,
    updateHierarchySpace,
)
from xmc import BayesianEstimator


class HierarchyOptimiser:
    """
    This class is in charge of deciding the best hierarchy. It is able to:
    - provide the best index set (according to its implement meaning of
            best) within adjustable constraints;
    - give the tolerance splitting associated to that choice (be it fixed or
            computed)
    - provide the sample numbers per Monte Carlo index that minimise the
            cost.

    It is assumed here that the hierarchy is optimised with respect to a single QoI.
    """

    indexSpace: List[int]
    """Lower and upper bounds on Monte Carlo indices in the hierarchy. Not yet implemented."""

    toleranceSplittingBounds: List[float]
    """Lower and upper bounds on the fraction of the total tolerance that the *statistical* error must satisfy. Set both bounds equal to keep this fraction constant."""

    _tolerances: List[float]
    """List of tolerances to be considered for optimisation, in order."""

    toleranceMultiplier: float
    """Number by which to multiply the current tolerance is a tolerance update is requested and there is no next tolerance. Must be in the interval ]0.0, 1.0]."""

    defaultHierarchy: List[List[Union[List[int], int]]]
    """Hierarchy chosen when no data is available to find the optimal one (e.g. at the first iteration)."""

    optimalIndexSet: Callable
    """Function used to find the optimal index set."""

    optimalSampleNumbers: Callable
    """Function used to find the optimal number of samples per Monte Carlo index."""

    updateHierarchySpace: Callable
    """Function used to update the space of acceptable hierarchies (e.g. bounds on indices and sample numbers.
    Not yet implemented."""

    toleranceSplittingMethod: Callable
    """Function used to optimise the splitting of the total tolerance between the different sources of error (e.g. statistical, bias)."""

    varianceBlender: BayesianEstimator
    """Object used to blend together model predictions and statistical estimates for index-wise data."""

    isVarianceBlended: bool
    """Whether model predictions and statistical estimates for index-wise data should be blended together."""

    def __init__(self, **keywordArgs):
        # Attributes
        self.indexSpace = keywordArgs.get("indexSpace", None)
        self.toleranceSplittingBounds = keywordArgs.get("toleranceSplittingBounds", None)
        self.defaultHierarchy = keywordArgs.get("defaultHierarchy", None)

        # Methods
        self.optimalIndexSet = dynamicImport(keywordArgs.get("optimalIndexSet"))
        self.optimalSampleNumbers = dynamicImport(keywordArgs.get("optimalSampleNumbers"))
        self.updateHierarchySpace = dynamicImport(
            keywordArgs.get("updateHierarchySpace", "xmc.tools.doNothing")
        )
        self.toleranceSplittingMethod = dynamicImport(
            keywordArgs.get("toleranceSplittingMethod", "xmc.tools.doNothing")
        )

        # Initialisation of mechanism to update tolerance
        self._tolerances = keywordArgs.get("tolerance", None)
        if isinstance(self._tolerances, float):
            self._tolerances = [self._tolerances]
        self.toleranceMultiplier = keywordArgs.get("toleranceMultiplier", 1.0)
        if self.toleranceMultiplier > 1.0 or self.toleranceMultiplier <= 0.0:
            raise ValueError(
                "toleranceMultiplier can only have a value in the interval ]0.0, 1.0]."
            )

        # TODO - Think of a better way to do this, improve naming
        self.varianceBlender = keywordArgs.get("varianceBlender", None)
        self.isVarianceBlended = keywordArgs.get("isVarianceBlended", False)

    # TODO This should be a static method
    def inputDictionaryTemplate(self):
        input_dict = {
            "oldHierarchy": None,
            "newSampleNumber": None,
            "errorParameters": None,
            "splittingParameter": None,
            "tolerances": None,
            "models": None,
            "parametersForModel": None,
            "estimations": None,
            "costModel": None,
            "costParameters": None,
            "costEstimations": None,
            "minimalSamplesPerIndex": 5,
            "blendedVariances": None,
            "defaultHierarchy": None,
        }
        return input_dict

    def toleranceSplitting(self, inputDict):
        # TODO THE FOLLOWING HOLDS GOOD ONLY FOR ONE TOLERANCE AND ONE SPLITTING PARAMETER
        # DO NOT USE if more are needed
        if len(inputDict["tolerances"]) > 1:
            raise ValueError(
                "toleranceSplitting() cannot handle more than one tolerance for now"
            )

        if self.toleranceSplittingBounds[0] == self.toleranceSplittingBounds[1]:
            inputDict["splittingParameter"] = self.toleranceSplittingBounds[0]
        else:
            splitting_parameter = self.toleranceSplittingMethod(inputDict)
            splitting_parameter = min(
                max(splitting_parameter, self.toleranceSplittingBounds[0]),
                self.toleranceSplittingBounds[1],
            )
            inputDict["splittingParameter"] = splitting_parameter
        return inputDict

    def updateTolerance(self):
        """
        Update current tolerance to next value in sequence
        """
        if not self._tolerances:
            pass

        if len(self._tolerances) > 1:
            del self._tolerances[0]
        else:
            self._tolerances[0] = self._tolerances[0] * self.toleranceMultiplier

    @property
    def tolerance(self) -> float:
        """Tolerance currently considered to optimise the hierarchy"""
        if self._tolerances:
            return self._tolerances[0]
        else:
            return None

    def optimalHierarchy(self, inputDict):
        """
        Compute the best hierarchy by individually computing the best index
        set and best number of samples
        """
        inputDict["tolerances"] = [self.tolerance]
        inputDict = self.toleranceSplitting(inputDict)
        new_ind = self.optimalIndexSet(inputDict)
        # Enforce index space bounds
        # 2nd condition necessary until MC indices are fixed
        if self.indexSpace and new_ind[0]:
            indices = [i for i in new_ind if i[0] <= max(self.indexSpace)]
        else:
            indices = new_ind
        # TODO - Think of a better way to do the following. Very fragile.
        if self.isVarianceBlended and inputDict["parametersForModel"][0]:
            inputDict["blendedVariances"] = self.varianceBlender.blend(indices, inputDict)
        sample_numbers = self.optimalSampleNumbers(inputDict, indices)
        sample_numbers = [
            max(inputDict["minimalSamplesPerIndex"], sample_numbers[i])
            for i in range(len(sample_numbers))
        ]

        old_hierarchy = inputDict["oldHierarchy"]
        new_hierarchy = []
        if len(old_hierarchy[0]) == 0:
            new_hierarchy = [indices, sample_numbers]
        else:
            new_hierarchy = mergeTwoListsIntoOne(indices, sample_numbers)

        return new_hierarchy
