from xmc.tools import dynamicImport
from xmc.methodDefs_estimationAssembler import assembleEstimation

# TODO Make this class a simple API definition, with only method and attribute declaration
class EstimationAssembler:
    """
    The class implements what is needed to produce the final estimations desired
    from the Monte Carlo method, given the statistical estimations provided by
    each index.
    """

    def __init__(self, **keywordArgs):
        # Attributes
        self.parameters = keywordArgs.get("parameters", None)

        # Methods
        self.assembleEstimation = dynamicImport(keywordArgs.get("assembleEstimation"))
