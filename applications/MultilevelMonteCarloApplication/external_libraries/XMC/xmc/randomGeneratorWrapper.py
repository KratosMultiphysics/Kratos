# Import Python libraries
import numpy.random as npr
from abc import ABCMeta, abstractmethod
from typing import List, Union, Sequence, Tuple

# XMC imports
from xmc.tools import dynamicImport
from xmc.methodDefs_randomGeneratorWrapper import generator

numpy_event_t = Union[float, List[float]]


class XMCRandomGenerator(metaclass=ABCMeta):
    """Abstract template for generators of random events.

    This is a template for classes responsible for the generation of random events.
    Every random generator class must inherit from it.
    Such a generator can be either standalone, or an interface to an external tool."""

    @abstractmethod
    def realisation(self):
        """Single realisation of a random event.

        Every child class _must_ write its own definition."""
        return


class RandomGeneratorWrapper(XMCRandomGenerator):
    """Wrapper around a function which generates a random event."""

    def __init__(self, **keywordArgs):
        self.parameters = keywordArgs["parameters"]
        self._generator = dynamicImport(keywordArgs["generator"])

    def realisation(self):
        return self._generator(*self.parameters)


class EventDatabase(XMCRandomGenerator):
    """Reads random events from a preset database.

    This class is derived from RandomGeneratorWrapper and generate events from a preset list,
    that is passed as argument to the class.
    """

    events: Sequence
    """Event 'database': sequence containing events in desired order."""

    _eventCounter: int
    """Internal counter to return events in order."""

    def __init__(self, events: Sequence, eventCounter: int = 0) -> None:
        """Create an EventDatabase object.

        :param events: Set directly to :py:attr:`~.events`.
        :param eventCounter: Set directly to :py:attr:`~._eventCounter`. Optional. Default: ``O``.
        """

        # Event database (iterable)
        self.events = events
        # Observe that self._eventCounter is re-initialized for each batch of the asynchronous algorithm.
        self._eventCounter = eventCounter

    def realisation(self) -> list:
        # Take modulo to loop over event collection
        position = self._eventCounter % len(self.events)
        # Get event
        event = self.events[position]
        # Increment counter
        self._eventCounter += 1
        # Return event
        return event


class NumPyGenerator(XMCRandomGenerator):
    """Interface between XMC and the random generator of NumPy."""

    rng: npr.Generator
    """Random number generator used to draw events."""

    distributions: List[Tuple[str, Union[list, tuple]]]
    """List of couples ``(method, arguments)``,
    where ``method`` is the name of a method of :py:attr:`~.rng` (i.e. of an object of type :py:class:`numpy.random.Generator`)
    and ``arguments`` is a sequence of positional arguments to be passed to that method.
    This information is used to generate random events.
    Example: ``[('normal', [0,1]), ('uniform', (1, 10, 3))]``."""

    def __init__(
        self, distributions: List[Sequence[Union[str, Sequence]]], seed: int = None
    ) -> None:
        """Constructor.

        :param distributions: Set to :py:attr:`~.distributions`. Sequence types will be converted to the ones expected for :py:attr:`~.distributions` if necessary; e.g. ``[('normal', [0,1]), ['uniform', [1, 10, 3]]]`` is acceptable and will be converted to ``[('normal', (0,1)), ('uniform', (1, 10, 3))]``

        :param seed: Used to set the initial state of the random generator. Optional, default to ``None`` (i.e. NumPy's own initialisation).
        :type seed: any type accepted by :py:meth:`numpy.random.default_rng`: ``None``, ``int``, sequences of ``int`` objects and various NumPy types (see the documentation of :py:meth:`~numpy.random.default_rng`).

        **Example.**

        Consider the following code::

            from xmc.randomGeneratorWrapper import NumPyGenerator
            rng = NumPyGenerator(seed=1234,
                distributions = [("normal", (0, 1)), ("uniform", (0, 1, 3))])

        The code above creates a :py:class:`~.NumPyGenerator` instance ``rng`` with seed ``1234``.
        The output of ``rng.realisation()`` will be a list of length 2, a single realisation of :math:`(n,u)\sim (\mathcal{N}(0,1), \mathcal{U}([0,1])^3)`.
        The random variable :math:`n` follows a normal law with mean 0 and standard deviation 1;
        the random variable :math:`u := (u_1, u_2, u_3)` has three independent components, each follows a uniform law over the interval :math:`[0,1]` i.e. :math:`\\forall i \\in\{1,2,3\}`, :math:`u_i \sim \mathcal{U}([0,1])`.
        """
        # Create random generator
        self.rng = npr.default_rng(seed)

        # Get distributions, that are names of methods (mandatory) and parameters
        self.distributions = []
        for dist in distributions:
            # Detect possible mistakes in format
            assert len(dist) == 2, (
                "Expected elements of 'distributions' to have length 2, " f"not {len(dist)}"
            )
            # dist[0] is the name of the method (str)
            # dist[1] is the sequence of argument, which we convert to tuple.
            self.distributions.append((dist[0], tuple(dist[1])))

    def realisation(self) -> List[numpy_event_t]:
        """Single realisation of random event.

        Calls all the distributions named in :py:attr:`~.distributions` of the NumPy random generator stored in :py:attr:`~.rng`."""
        return [self.generatorCall(method, args) for method, args in self.distributions]

    def generatorCall(self, methodName: str, arguments: tuple) -> numpy_event_t:
        """Single call to a method of the random generator.

        :param methodName: name of the method to call. Must be a callable attribute of :py:attr:`~.rng`, i.e. of an object of class :py:class:`numpy.random.Generator`.
        :param arguments: sequence of arguments to be passed to the callable attribute named by ``methodName``.
        :return: the random event generated by this call, converted to a built-in Python type beforehand.
        """
        rng_method = getattr(self.rng, methodName)
        event = rng_method(*arguments)
        # Convert to Python built-in types
        if not isinstance(event, float):
            event = event.astype(float).tolist()
        return event
