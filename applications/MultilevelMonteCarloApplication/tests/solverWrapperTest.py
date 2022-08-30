"""Test cases for SolverWrapper classes."""

# Import python class test
import unittest

# External tools
import numpy as np
import json
from typing import Callable, Type, Optional, List, Tuple

# Import XMC class
from xmc.solverWrapper import SolverWrapper
from xmc.classDefs_solverWrapper.oscillator import VanDerPol


class TestSolver:
    """Abstract test case templates for solver classes"""

    class Base(unittest.TestCase):
        """Template test case for SolverWrapper instances.

        To be inherited. When creating a sub-class, overload the class method ``setUpClass`` to set the attribute ``solver`` to the class to be tested. See  :py:class:`~.TestSolverWrapper` for an example.
        """

        solverClass: Type[SolverWrapper]
        """Class tested"""

        classDefaults: dict
        """Dictionary whose keys are all *optional* keyword arguments passed to the class constructor to set the attribute of *same name* as the keyword. The value of each entry is the default value of this attribute. Example: class ``MySolver`` accepts keyword argument ``thisAttribute`` to set the attribute ``MySolver.thisAttribute``. This keyword argument is optional, and the default value of ``MySolver.thisAttribute`` is ``that_value``. Then ``classDefaults[thisAttribute] = that_value``."""

        testDefaults: dict
        """Dictionary whose keys are all *mandatory* keyword arguments passed to the class constructor. The union of :py:attr:`~.testDefaults` and :py:attr:`~.classDefaults` must be enough to create an instance of :py:attr:`~.solverClass`."""

        solveArgs: list
        """List of positional arguments that can be successfully passed to ``solverClass.solve()``; i.e. ``self.solverClass.solve(*self.solveArgs)`` must not fail."""

        @classmethod
        def setUpClass(cls) -> None:
            """Sets solver class and default values"""

            # Set to callable placeholder here
            cls.solverClass = lambda *args, **kwargs: None

            # Default values of the class
            cls.classDefaults = {}
            # Default values used in this test case
            # for mandatory arguments without default values
            cls.testDefaults = {"solverWrapperIndex": (1,), "outputDimension": 1}
            cls.solveArgs = [1]

        def solver(self, *args, **kwArgs) -> SolverWrapper:
            """Instance of the solver to be tested.

            Parameters are passed to the class constructor.

            :return: solver object of tested class
            :rtype: :py:attr:`~.solverClass`
            """
            # Sets necessary arguments to test default values, if not given
            for key, val in self.testDefaults.items():
                # TODO Poor workaround until this attribute named is changed
                if key == "solverWrapperIndex":
                    kwArgs["index"] = kwArgs.get(key, val)
                else:
                    kwArgs[key] = kwArgs.get(key, val)
            return self.solverClass(*args, **kwArgs)

        def test_init(self) -> None:
            """Tests the creation of a solver instance"""
            with self.subTest(msg="Default values: creation without user input"):
                s = self.solver()
                for attr, val in dict(self.classDefaults, **self.testDefaults).items():
                    with self.subTest(attribute=attr, expected_value=val):
                        self.assertTrue(hasattr(s, attr))
                        if val is None:
                            self.assertIsNone(getattr(s, attr), val)
                        else:
                            self.assertEqual(getattr(s, attr), val)

        def test_solve(self) -> None:
            """Tests the solve method"""
            s = self.solver()
            with self.subTest(msg="Test solve returns correct types"):
                sol, cost = s.solve(*self.solveArgs)
                # Check cost
                self.assertIsInstance(cost, float)
                self.assertGreater(cost, 0)
                # Check solution
                self.assertIsInstance(sol, list)
                if isinstance(s.outputDimension, int):
                    self.assertEqual(len(sol), s.outputDimension)
                else:
                    self.assertEqual([len(e) for e in sol], s.outputDimension)


class SolverWrapperTest(TestSolver.Base):
    """Test case for the SolverWrapper class itself.

    Identical to the template :py:class:`.TestSolver.Base`. To be run as a unit test case.
    """

    @classmethod
    def setUpClass(cls) -> None:
        """Sets solver class."""

        super().setUpClass()
        cls.solverClass = SolverWrapper

        # Skip irrelevant tests
        cls.test_solve = unittest.skip("Test irrelevant for template class.")(cls.test_solve)


class VanDerPolTest(TestSolver.Base):
    """Test case for the VanDerPol solver class.

    Based on the template :py:class:`.TestSolver.Base` for SolverWrapper instances. To be run as a unit test case.
    """

    referenceTrajectoryFile: str
    """Relative path to JSON file containing reference solution. Used to test resolution methods."""

    @classmethod
    def setUpClass(cls) -> None:
        """Sets solver class and defaults.

        Overloads :py:meth:`.TestSolver.Base.setUpClass`.
        """

        super().setUpClass()
        cls.solverClass = VanDerPol
        del cls.testDefaults["outputDimension"]
        del cls.testDefaults["solverWrapperIndex"]
        cls.testDefaults["timestep"] = 1.0
        cls.testDefaults["duration"] = 100.0
        cls.classDefaults["initial"] = (1.0, 1.0)
        cls.classDefaults["damping"] = 0.0
        cls.classDefaults["forcingAmplitude"] = 1.0
        cls.classDefaults["offset"] = 0.0
        nb_steps = int(cls.testDefaults["duration"] / cls.testDefaults["timestep"])
        cls.solveArgs = [[1 for _ in range(0, nb_steps)]]
        cls.referenceTrajectoryFile = "parameters/vanderpol_reference_solution.json"

    def test_init(self) -> None:
        """Tests the creation of a VanDerPol instance.

        Overloads :py:meth:`.TestSolver.Base.test_init`.
        """
        # Overload
        super().test_init()

        solution_values = (
            "average",
            "average_squared",
            "variance",
            "deviation",
            np.linalg.norm,
        )
        trajectory = [1.0, -2.0, 3.0, 4.0]
        for sol in solution_values:
            with self.subTest(msg="Non-defaults solution choice", solution=sol):
                vdp = self.solver(solution=sol)
                # What follows goes beyond testing just init
                # but helps ensuring that a valid object was assigned to vdp.processTrajectory
                s = vdp.processTrajectory(trajectory)
                self.assertIsInstance(float(s), float)

        # Test refinement of timestep
        r = 3
        i = 2
        with self.subTest(
            msg="Refinement per solver index", timestep_multiplier=r, solver_index=i
        ):
            vdp_0 = self.solver(index=(0,))
            vdp_1 = self.solver(index=(i,), timestepDivider=r)
            mesh_0 = vdp_0.mesh
            mesh_1 = vdp_1.mesh
            self.assertEqual(len(mesh_0) - 1, (len(mesh_1) - 1) / r ** i)
            self.assertEqual(mesh_0[0], mesh_1[0])
            self.assertEqual(mesh_0[-1], mesh_1[-1])
        r = 3.0
        with self.subTest(
            msg="Refinement per solver index", timestepMultiplier=r, solver_index=i
        ):
            vdp = self.solver(index=(i,), timestepDivider=r)
        r = 3.5
        with self.subTest(
            msg="Refinement per solver index", timestepMultiplier=r, solver_index=i
        ):
            with self.assertWarns(UserWarning):
                # Must return a warning because timestepDivider is not an integer
                vdp = self.solver(index=(i,), timestepDivider=r)
            expected_timestep = self.testDefaults["timestep"] / int(r) ** i
            self.assertEqual(vdp.timestep, expected_timestep)

    def test_harmonic_trajectory(self) -> None:
        """Tests VanDerPol trajectory method in harmonic case.

        The reference is the analytical solution of a free harmonic oscillator. This correspond to a Van Der Pol oscillator with no damping and no forcing.
        The VanDerPol solver uses the Euler scheme, which does not conserve energy. As a result, the error between the computed trajectory and the reference grows with time. It is kept reasonably small by a short duration and small timesteps.
        """
        # Non-trivial initial conditions
        initial = (1.5, 0.8)
        # Create solver
        # The longer the duration, the smaller the timestep
        vdp = self.solver(
            duration=4 * np.pi,
            timestep=0.001,
            damping=0.0,
            forcingAmplitude=0.0,
            initial=initial,
        )
        # Compute trajectory
        trajectory = vdp.trajectory()
        # Compute reference trajectory
        times = vdp.mesh
        reference = VanDerPolTest.harmonic_trajectory(times, initial, 1.0)
        # Error between solver result and reference
        error = np.linalg.norm(trajectory - reference, ord=np.inf)
        # Allow  large tolerance because Euler scheme does not conserve energy
        self.assertAlmostEqual(error, 0.0, delta=0.01)

    def test_reference_solution(self) -> None:
        """Tests VanDerPol solution against saved reference.

        The reference trajectory was generated by the same code through :py:meth:`~.VanDerPolTest.generate_reference_solution`. This is not a validation of the solver, but a mean to detect regressions.
        """
        # import reference data
        # The data is assumed to be small, so we store it all
        with open(self.referenceTrajectoryFile, "r") as import_file:
            reference = json.load(import_file)
        vdp = self.solver(**reference["solverArguments"])
        trajectory = vdp.trajectory(reference["forcing"])
        trajectory_error = np.linalg.norm(
            trajectory - np.array(reference["trajectory"]), ord=2
        )
        self.assertAlmostEqual(trajectory_error, 0.0, delta=10 ** -11)
        solution = vdp.solve(reference["forcing"])[0][0]
        solution_error = np.abs(solution - np.array(reference["solution"]))
        self.assertAlmostEqual(solution_error, 0.0, delta=10 ** -11)

    @staticmethod
    def harmonic_trajectory(
        times: List[float],
        initial: Optional[Tuple[float, float]] = (1.0, 0.0),
        pulsatance: Optional[float] = 1.0,
    ) -> List[float]:
        """Trajectory of a simple harmonic oscillator.

        Amplitude and phase are deduced from given initial conditions.

        :param times: temporal mesh, as from :py:attr:`VanDerPol.mesh`
        :param initial: initial position and speed, as in :py:attr:`VanDerPol.initial`. Default: ``(1.0,0.0)``
        :param pulsatance: angular frequency of the harmonic oscillator, i.e. :math:`\\text{pulsatance} = 2 \\pi \\times \\text{frequency}`. Default: ``1.0``.
        :return: positions of the oscillator at times ``times``, as from :py:meth:`VanDerPol.trajectory`.
        """
        if initial[0] == 0:
            phase = np.pi / 2
            amplitude = -initial[1] / (pulsatance * np.sin(phase))
        else:
            phase = np.arctan(-initial[1] / (pulsatance * initial[0]))
            amplitude = initial[0] / np.cos(phase)
        return amplitude * np.cos(pulsatance * np.array(times) + phase)

    @staticmethod
    def generate_reference_solution(
        filePath: str, solver: VanDerPol, forcing: Optional[List[float]] = None
    ) -> None:
        """Produce reference solution data required by :py:meth:`~.VanDerPolTest.test_reference_solution`.

        The reference solution data is computed using a :py:class:`VanDerPol` instance. Both the solution data and the parameters to reproduce it are written to a JSON file.

        :param filePath: path to the file to write. File will be completey overwritten if it exists, and created if it does not.
        :param solver: :py:class:`VanDerPol` instance to use. Its :py:attr:`~VanDerPol.processTrajectory` attribute will be overwritten to ``np.mean``, because it cannot be serialised to be stored in the JSON file.
        :param forcing: event to use as input to :py:meth:`VanDerPol.solve`. Default: random generation from a standard normal distribution.
        """
        if forcing is None:
            forcing = np.random.default_rng().normal(0.0, 1.0, size=len(solver.mesh)).tolist()
        # Override
        solver.processTrajectory = np.mean
        # Get trajectory
        trajectory = solver.trajectory(forcing)
        # Get time average
        # The output of solve is ([time_average],cost)
        solution = solver.solve(forcing)[0][0]
        # Get dictionary
        solver_dict = solver.__dict__
        # Replace (deprecated) attribute solverWrapperIndex by corresponding input keyword
        solver_dict["index"] = solver_dict.pop("solverWrapperIndex")
        # Remove unserialisable attributes
        del solver_dict["processTrajectory"]
        # Assemble data to be exported
        export = {
            "solverArguments": solver_dict,
            "forcing": forcing,
            "trajectory": trajectory,
            "solution": solution,
        }
        # Write to file
        with open(filePath, "w+") as export_file:
            json.dump(export, export_file)
