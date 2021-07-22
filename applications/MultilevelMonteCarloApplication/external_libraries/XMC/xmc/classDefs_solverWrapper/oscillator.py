"""Solvers for oscillator problems."""

import time
import numpy as np
from warnings import warn

from xmc.solverWrapper import SolverWrapper

# Types
from typing import Optional, Callable, List, Tuple

# TODO: Annotations below are not displayed by Sphinx (autodoc)

event_t = List[float]
"""Type of a random event expected by ``solve``."""

output_t = List[float]
"""Type of a quantity of interest returned by the solver."""

cost_t = float
"""Type of the cost of a resolution, returned by the solver."""


class VanDerPol(SolverWrapper):
    """Python solver for the Van der Pol oscillator problem with random forcing."""

    timestep: float
    """size of step between two points in the uniform temporal mesh."""

    duration: float
    """total duration of the simulation."""

    initial: Tuple[float, float]
    """initial values of position and speed of the oscillator."""

    damping: float
    """damping parameter of the oscillator."""

    forcinAmplitude: float
    """Multiplicator factor to be applied to the forcing."""

    offset: float
    """Constant value to be substracted to the trajectory *before* applying :py:attr:`~.processTrajectory`."""

    processTrajectory: Callable
    """Function to be applied to oscillator trajectory to get desired solution.

    :param list[float] trajectory: oscillator trajectory as returned by :py:meth:`~.trajectory`.
    :return: desired solution, as float (or an object which can be converted by ``float()``).
    :rtype: float
    """

    def __init__(self, **kwArgs) -> None:
        """Constructor.

        :keyword float timestep: Used to set :py:attr:`~.timestep`, along with ``timestepDivider`` and ``index``.

        :keyword timestepDivider: Factor used to change the timestep according to the index of the solver. :py:math:`\\text{self.timestep} = \\text{timestep} \\times \\text{timestepDivider} ^ {-\\text{index}}`. Should be an ``int``, but integer floats are accepted. Default: ``2``.
        :type timestepDivider: int or None

        :keyword float duration: Sets :py:attr:`~.duration`.

        :keyword initial: Sets :py:attr:`~.initial`. Default: ``(1.0, 1.0)``
        :type initial: tuple(float,float) or None

        :keyword damping:  Sets :py:attr:`~.damping`. Default: ``0.0``
        :type damping: float or None

        :keyword forcingAmplitude: Sets :py:attr:`~.forcingAmplitude`. Default: ``1.0``
        :type forcingAmplitude: float or None

        :keyword offset:  Sets :py:attr:`~.offset`. Default: ``0.0``
        :type offset: float or None

        :keyword solution: Used to set :py:attr:`~.processTrajectory`. If it is a callable object, it is set directly to :py:attr:`~.processTrajectory` and must therefore conform to its specifications. If it is a string, it must be one of the preset values: ``"average"``, ``"average_squared"``, ``"variance"``, ``"deviation"``. Default: ``"average"``
        :type solution: str or callable

        :keyword index: Sets :py:attr:`~.solverWrapperIndex`. See :py:attr:`SolverWrapper.solverWrapperIndex <xmc.solverWrapper.SolverWrapper.solverWrapperIndex>`. Default: ``(0,)``
        :type index: tuple[int] or None

        :keyword outputDimension: Sets :py:attr:`~.outputDimension`. See :py:attr:`SolverWrapper.outputDimension <xmc.solverWrapper.SolverWrapper.outputDimension>`. Default: ``1``
        :type outputDimension: int or list[int] or None
        """
        # Provide default values to make parameters optional
        kwArgs["index"] = kwArgs.get("index", (0,))
        kwArgs["outputDimension"] = kwArgs.get("outputDimension", 1)
        super().__init__(**kwArgs)

        # Refine timestep according to index
        base_timestep = kwArgs.get("timestep")
        r = kwArgs.get("timestepDivider", 2)
        # Ensure that the multiplier is an integer (even if it is of type float)
        if isinstance(r, float) and r != int(r):
            warn("timestepDivider must be an integer. I will truncate it", UserWarning)
            r = int(r)
        self.timestep = base_timestep * r ** -self.solverWrapperIndex[0]  # dt

        # Other problem parameters (deterministic)
        self.duration = kwArgs.get("duration")  # T
        self.initial = kwArgs.get("initial", (1.0, 1.0))  # x0, dx0
        self.damping = kwArgs.get("damping", 0.0)  # mu
        self.forcingAmplitude = kwArgs.get("forcingAmplitude", 1.0)  # sigma or tau
        self.offset = kwArgs.get("offset", 0.0)

        # Process choice of output solution
        output_choice = kwArgs.get("solution", "average")
        if callable(output_choice):
            self.processTrajectory = output_choice
        elif output_choice == "average":
            self.processTrajectory = np.mean
        elif output_choice == "average_squared":
            self.processTrajectory = VanDerPol._trajectoryAverageSquared
        elif output_choice == "variance":
            self.processTrajectory = np.var
        elif output_choice == "deviation":
            self.processTrajectory = np.std
        else:
            raise ValueError(
                f"Unexpected value for keyword argument solution: {output_choice}"
            )

    def solve(self, event: event_t) -> Tuple[output_t, cost_t]:
        """Solves the oscillator problem for the given event and returns the desired quantity.

        The trajectory is computed by the :py:meth:`~.VanDerPol.trajectory` method.

        :param event: External forcing. If it is finer than the temporal mesh, it will be :py:meth:`coarsened <.VanDerPol.coarsenResolution>`.
        :return: The temporal average of the trajectory (in a list) and the time taken by the method.
        """
        # Start counting time
        start_time = time.time()

        # Compute trajectory
        random_forcing = self.coarsenResolution(event)
        trajectory = self.trajectory(random_forcing)
        # Offset trajectory
        trajectory = (np.array(trajectory) - self.offset).tolist()
        # Compute requested solution, stored as type output_t
        value = [float(self.processTrajectory(trajectory))]

        # Measure computational time (real time)
        cost = time.time() - start_time
        return value, cost

    def trajectory(self, forcing: Optional[List[float]] = None) -> List[float]:
        """Trajectory of the oscillator, with given forcing.

        Solves the oscillator problem with forward Euler scheme and uniform temporal mesh.

        :param forcing: list of values of external forcing, as discretised on :py:attr:`~.VanDerPol.mesh` (must have same length). If `None`, no forcing is applied.
        :return: values of position at times of :py:attr:`~.VanDerPol.mesh`.
        """
        # Shorthands
        dt = self.timestep
        nb_dt = self._numberTimesteps()
        damping = self.damping
        # Process optional input argument
        if forcing:
            forcing_factor = self.forcingAmplitude * np.sqrt(dt)
        else:
            forcing_factor = 0
            forcing = np.zeros(nb_dt)
        x = np.zeros(nb_dt)
        dx = np.zeros(len(x))
        # Initialize
        x[0], dx[0] = self.initial

        # Solve iteratively
        for n in range(nb_dt - 1):
            x[n + 1] = x[n] + dt * dx[n]
            dx[n + 1] = (
                dx[n]
                + dt * (damping * dx[n] * (1 - x[n] ** 2) - x[n])
                + forcing_factor * forcing[n]
            )

        return x.tolist()

    def coarsenResolution(self, timeValues: event_t) -> event_t:
        """Reduce temporal resolution of a discrete time process.

        It assumes that the array received corresponds to a temporal mesh which is a uniform refinement of that of this object. The coarsening is achieved by summation of intermediate values (thus treated as increments).
        This is intended for random forcing received as input to :py:meth:`~.VanDerPol.solve`.

        :param timeValues: values of length greater than that of :py:attr:`~.VanDerPol.mesh`. If ``len(timeValues) <= len(self.mesh)-1``, ``timeValues is returned unchanged``. The comparison is done with ``len(self.mesh)-1`` and not ``len(self.mesh)`` because :py:attr:`~.VanDerPol.mesh` includes the initial timestep, while ``timeValues`` must not.
        :return: list of values matching the length of :py:attr:`~.VanDerPol.mesh`.
        """
        # Check if coarsening is needed
        nb_steps = self._numberTimesteps() - 1  # exclude initial timestep
        if nb_steps >= len(timeValues):
            return timeValues

        # Compute ratio of resolutions in temporal domain
        ratio = int(np.floor(len(timeValues) / nb_steps))
        # Coarsen resolution
        r = np.sqrt(ratio)
        coarsenedValues = [
            sum(timeValues[i : i + ratio]) / r for i in range(0, len(timeValues), ratio)
        ]

        return coarsenedValues

    def _numberTimesteps(self) -> int:
        """Number of timesteps used in the resolution (including initial)."""
        return int(np.ceil(self.duration / self.timestep)) + 1

    @property
    def mesh(self) -> List[float]:
        """The temporal mesh used.

        The mesh uniform and defined by attributes :py:attr:`~.VanDerPol.timestep` and :py:attr:`~.VanDerPol.duration`. Change these to change the mesh.
        """
        # NB arange exclude upper bound, so timestep must be added to second argument
        return np.arange(0, self.timestep * self._numberTimesteps(), self.timestep).tolist()

    # Solution methods
    @staticmethod
    def _trajectoryAverageSquared(trajectory: List[float]) -> float:
        """Average of the square of given trajectory.

        :param trajectory: oscillator trajectory as returned by :py:meth:`~.trajectory`.
        """
        trajectory = np.array(trajectory) ** 2
        return float(np.mean(trajectory))
