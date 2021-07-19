# XMC imports
from xmc.tools import dynamicImport
from xmc.methodDefs_randomGeneratorWrapper import generator

# TODO RandomGeneratorWrapper must have a 'seed' attribute for reproducibility. Its generator method must support it.


class RandomGeneratorWrapper:
    """
    This is a class is responsible for the generation of random inputs: scalars,
    fields, whatever. It is very similar to SolverWrapper: it handles calls to
    whatever tool (most likely external) is chosen for this task. Therefore, this
    is a barebone template from which tool-specific implementation will be
    fashioned, either by inheritance or pointers.
    """

    def __init__(self, **keywordArgs):
        self.parameters = keywordArgs.get("parameters")
        self._generator = dynamicImport(keywordArgs.get("generator"))

    def realisation(self):
        return self._generator(*self.parameters)

class EventDatabase(RandomGeneratorWrapper):
    """
    This class is derived from RandomGeneratorWrapper and generate events from a preset list,
    that is passed as argument to the class.
    """

    def __init__(self, **keywordArgs) -> None:
        """
        Create an EventDatabase object.

        Keyword arguments
        ---------------

        events: list.
            Iterable containing events in desired order.
        """

        # Provide safe default values
        keywordArgs["parameters"] = []
        keywordArgs["generator"] = "xmc.tools.doNothing"
        # Parent constructor
        super().__init__(**keywordArgs)
        # Event database (iterable)
        self.events = keywordArgs.get("events")
        # Internal counter to return events in order
        # Observe that self._eventCounter is re-initialized for each batch of the asynchronous algorithm.
        self._eventCounter = 0

    def realisation(self) -> list:
        # Take modulo to loop over event collection
        position = self._eventCounter % len(self.events)
        # Get event
        event = self.events[position]
        # Increment counter
        self._eventCounter += 1
        # Return event
        return event