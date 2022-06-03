# Import Python libraries
import numpy as np

# XMC imports
from xmc.tools import dynamicImport
from xmc.methodDefs_multiCriterion import flag, interpreter


class MultiCriterion:
    """
    This class handles any number of criteria, the associated tolerances (with
    possible splitting), and the combination of these criteria into an interpred output.
    """

    def __init__(self, **keywordArgs):
        # Attributes
        self.criteria = keywordArgs.get("criteria")
        self.inputsForCriterion = keywordArgs.get("inputsForCriterion")
        self._interpreter = dynamicImport(keywordArgs.get("interpreter"))
        self.splitCriteria = keywordArgs.get("splitCriteria", None)
        self.toleranceToSplit = keywordArgs.get("toleranceToSplit", None)
        self.flagStructure = interpreter.interpretationStructure

        # Methods
        self._flagFunction = dynamicImport(
            keywordArgs.get("flag", "xmc.methodDefs_multiCriterion.flag.plainFlag")
        )

    def __getstate__(self):
        # Captures what is normally pickled
        state = self.__dict__.copy()
        # Replace the PyCOMPSs-decorated entries with explicit label
        for attribute, value in state.items():
            if hasattr(value, "__module__") and "pycompss" in value.__module__:
                state[attribute] = "unassigned_task_decorator"
        # what we return here will be stored in the pickle
        return state

    # TODO Find a solution to re-assign original PyCOMPS-decorated attributes
    # or at least an undecorated version.
    # See reference below for a simple way to do that (but it is intrusive to PyCOMPSs)
    # https://stackoverflow.com/a/33024739
    # This is not currently necessary, just more robust.
    #    def __setstate__(self,newState):
    # Re-create desired instance
    # ...
    # re-instate our __dict__ state from the pickled state
    #        self.__dict__.update(newState)

    def tolerances(self, criteriaReferences=None):
        """
        Returns the tolerances of the requested criteria.
        """
        # TODO Is it possible to set default as criteria_references=range(0,len(self.criteria)-1) instead of using this conditional structure?
        if criteriaReferences is None:
            criteriaReferences = range(len(self.criteria))
        tolerance_values = []
        for i in criteriaReferences:
            # TODO Do not add None values to tolerance_values
            # It is perfectly normal that some MonoCriterion objects have tolerance None.
            tolerance_values.append(self.criteria[i].tolerance)
        return tolerance_values

    def splittingParameter(self):
        """
        Returns the splitting parameter currently applied
        """
        return self.criteria[self.splitCriteria[0]].tolerance / self.toleranceToSplit

    def setTolerance(self, criterionReferences, toleranceValues):
        for i in range(len(criterionReferences)):
            self.criteria[criterionReferences[i]].tolerance = toleranceValues[i]

    def splitTolerance(self, splittingParameter):
        if self.toleranceToSplit is not None:
            split_tolerances = self.toleranceToSplit * np.array(
                [splittingParameter, 1 - splittingParameter]
            )
            self.setTolerance(self.splitCriteria, split_tolerances)
        else:
            pass

    def updateTolerance(self, criteriaToUpdate=None):
        """
        Update the tolerances of self.criteria entries specified by criteriaToUpdate
        """
        if criteriaToUpdate is None:
            criteriaToUpdate = range(len(self.criteria))

        for coord in criteriaToUpdate:
            if len(self.criteria[coord].tolerances) == 0:
                ValueError(
                    "stoppingCriterion.criteria[", coord, "] has no tolerance to update"
                )
            elif len(self.criteria[coord].tolerances) == 1:
                self.criteria[coord].tolerance = self.criteria[coord].tolerances[0]
            else:
                self.criteria[coord].tolerance = self.criteria[coord].tolerances.pop(0)

    def flag(self, values):
        """
        Return the output of the criterion.
        It takes the expected values as input, evaluate each elementary criterion on them and combine their boolean outputs into a dictionary flag.
        This is currently a wrapper for the protected _flagFunction attribute.
        """
        return self._flagFunction(
            values, self.criteria, self.inputsForCriterion, self._interpreter
        )
