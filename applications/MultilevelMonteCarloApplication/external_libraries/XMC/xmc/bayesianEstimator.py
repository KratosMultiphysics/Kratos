# XMC imports
import xmc.statisticalEstimator as se
from xmc.tools import dynamicImport
from xmc.methodDefs_bayesianEstimator import blend


class BayesianEstimator(se.StatisticalEstimator):
    """
    This class is for Bayesian estimation, where an a priori estimation is updated
    by observations.
    """

    def __init__(self, **keywordArgs):
        # Parent constructor
        super().__init__(**keywordArgs)

        # Attributes
        self.parameters = keywordArgs.get("parameters", None)

        # Methods
        self._blendFunction = dynamicImport(
            keywordArgs.get("blendFunction", "xmc.bayesianEstimator.blend.noBlending")
        )

    def blend(self, newLevels, inputDict):
        return self._blendFunction(self.parameters, newLevels, inputDict)
