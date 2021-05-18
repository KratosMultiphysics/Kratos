# XMC imports
from xmc.tools import dynamicImport
from xmc.methodDefs_modelEstimator import update, valueForParameters


class ModelEstimator:
    """
    Generic class for model with a known analytical expression, dependent on
    parameters to be fitted on data. A data point is an observed (antecedent ;
    image) couple.
    Some methods descriptions.
    - value: return model evaluation at given point.
    - valueForParameters: Returns the model evaluation at given point with given parameter values.
    - fitDistance: computes the distance between two sets of data points; this is the distance used to fit the model to new observations.
    """

    def __init__(self, **keywordArgs):
        # Attributes
        self.parameters = keywordArgs.get("parameters", None)
        self.oldParameters = keywordArgs.get("oldParameters", None)
        self.isDataStored = keywordArgs.get("isDataStored", False)
        self.data = keywordArgs.get("data", [])

        # Methods
        self.updater = dynamicImport(keywordArgs.get("updater"))
        self._valueForParameters = dynamicImport(keywordArgs.get("valueForParameters"))
        self._fitDistance = dynamicImport(keywordArgs.get("fitDistance", None))

    def value(self, antecedentList):
        return self._valueForParameters(self.parameters, antecedentList)

    def update(self, *data):
        self.oldParameters = self.parameters
        self.parameters = self.updater(*data)
