# XMC imports
from xmc.tools import dynamicImport

# TODO RandomGeneratorWrapper must have a 'seed' attribute for reproducibility. Its generator method must support it.

class RandomGeneratorWrapper():
    """
    This is a class is responsible for the generation of random inputs: scalars,
    fields, whatever. It is very similar to SolverWrapper: it handles calls to
    whatever tool (most likely external) is chosen for this task. Therefore, this
    is a barebone template from which tool-specific implementation will be
    fashioned, either by inheritance or pointers.
    """

    def __init__(self, **keywordArgs):
        self.parameters = keywordArgs.get('parameters')
        self._generator = dynamicImport(keywordArgs.get('generator'))

    def realisation(self):
        return self._generator(*self.parameters)
