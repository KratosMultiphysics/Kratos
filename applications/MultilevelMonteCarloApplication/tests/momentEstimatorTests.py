# Import python class test
import unittest

# Import python libraries
from json import load
from copy import deepcopy
from itertools import product as itproduct

# Math tools
from math import factorial
from xmc.tools import doubleFactorial, binomial, returnInput_Task, normalInverseCDF
import numpy as np

# Import xmc classes
from xmc.momentEstimator import (
    MomentEstimator,
    CombinedMomentEstimator,
    MultiMomentEstimator,
    MultiCombinedMomentEstimator,
)
import xmc.methodDefs_momentEstimator.computeCentralMoments as mdccm

# Import exaqute for parallel tests
from exaqute import (
    task,
    get_value_from_remote,
    Type,
    COLLECTION_IN,
    Depth,
)

# Define types
from typing import List, DefaultDict, Union

PowerSumsDict = DefaultDict[str, List[float]]
SampleArray = List[List[List[float]]]
HStatistics = Union[float, List[float]]


class TestMomentEstimator(unittest.TestCase):
    def test_update(self):
        list_values = [[1.0], [2.0], [3.0], [4.0], [5.0], [6.0]]
        true_power_sums = [
            [21.0],
            [91.0],
            [441.0],
            [2275.0],
            [12201.0],
            [67171.0],
            [376761.0],
            [2142595.0],
            [12313161.0],
            [71340451.0],
        ]

        # read parameters
        with open("parameters/parameters_test_momentEstimator.json", "r") as parameter_file:
            parameters = load(parameter_file)

        for order in [1, 2, 5]:
            parameters["momentEstimatorInputDict"]["order"] = order
            parameters["momentEstimatorInputDict"]["updatedPowerSums"] = (
                "xmc.methodDefs_momentEstimator.updatePowerSums.updatePowerSumsOrder"
                + str(2 * order)
                + "Dimension0"
            )  # required order is 2 * order

            # build momentEstimator class
            test_me = MomentEstimator(**parameters["momentEstimatorInputDict"])

            # update power sums
            for value in list_values:
                test_me.update([value])

            # test update sample number
            self.assertEqual(test_me._sampleCounter, len(list_values))

            # test update power sums
            for i in range(2 * test_me.order):
                self.assertEqual(test_me.powerSums[i], true_power_sums[i])

    def test_value(self):
        list_values = [[1.0], [2.0], [3.0], [4.0], [5.0]]
        true_estimation_order_1 = 3.0
        true_estimation_order_2 = 2.5
        true_error_order_1 = 0.5
        true_error_order_2 = 0.975

        # read parametrs
        with open("parameters/parameters_test_momentEstimator.json", "r") as parameter_file:
            parameters = load(parameter_file)

        # build momentEstimator class
        test_me = MomentEstimator(**parameters["momentEstimatorInputDict"])

        # update power sums
        test_me.update(list_values)

        # test
        xmc_estimation_order_1 = test_me.value(
            order=1, isCentral=True, isErrorEstimationRequested=False
        )
        xmc_estimation_order_2 = test_me.value(
            order=2, isCentral=True, isErrorEstimationRequested=False
        )
        xmc_error_order_1 = test_me.value(
            order=1, isCentral=True, isErrorEstimationRequested=True
        )
        xmc_error_order_2 = test_me.value(
            order=2, isCentral=True, isErrorEstimationRequested=True
        )
        self.assertEqual(xmc_estimation_order_1, true_estimation_order_1)
        self.assertEqual(xmc_estimation_order_2, true_estimation_order_2)
        self.assertEqual(xmc_error_order_1, true_error_order_1)
        self.assertAlmostEqual(
            xmc_error_order_2, true_error_order_2
        )  # fail if the two objects are unequal
        # as determined by their difference rounded
        # to the given number of decimal places(default 7)
        # (according to unittest docs)


class TestCombinedMomentEstimator(unittest.TestCase):
    def test_update(self):
        Q = 2.0
        list_values = [
            [
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                10,
            ],
            [
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                10,
            ],
            [
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                10,
            ],
            [
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                10,
            ],
            [
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                10,
            ],
            [
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                [10 * Q],
                [10 * Q * Q],
                10,
            ],
        ]

        true_power_sums = [
            [120.0],
            [240.0],
            [120.0],
            [240.0],
            [120.0],
            [240.0],
            [120.0],
            [240.0],
            [120.0],
            [240.0],
        ]

        # read parameters
        with open(
            "parameters/parameters_test_combinedMomentEstimator.json", "r"
        ) as parameter_file:
            parameters = load(parameter_file)

        for order in [1, 5]:
            parameters["momentEstimatorInputDict"]["order"] = order
            # build momentEstimator class
            test_me = CombinedMomentEstimator(**parameters["momentEstimatorInputDict"])

            # update power sums
            for value in list_values:
                test_me.update([[value]])

            # test update sample number
            self.assertEqual(test_me._sampleCounter, 60)

            # test update power sums
            for i in range(2 * test_me.order):
                self.assertEqual(test_me.powerSums[i], true_power_sums[i])

            S1 = test_me.powerSums[0][0]
            S2 = test_me.powerSums[1][0]
            h1 = mdccm.computeCentralMomentsOrderOneDimensionZero(S1, test_me._sampleCounter)
            h2 = mdccm.computeCentralMomentsOrderTwoDimensionZeroBiased(
                S1, S2, test_me._sampleCounter
            )

            self.assertEqual(h1, 2.0)
            self.assertEqual(h2, 0.0)

    def test_value(self):
        Q = 2.0
        list_values = [
            [[10 * Q], [10 * Q * Q], 10],
            [[10 * Q], [10 * Q * Q], 10],
            [[10 * Q], [10 * Q * Q], 10],
            [[10 * Q], [10 * Q * Q], 10],
            [[10 * Q], [10 * Q * Q], 10],
            [[10 * Q], [10 * Q * Q], 10],
        ]

        true_estimation_order_1 = 2.0
        true_estimation_order_2 = 0.0
        true_error_order_1 = 0.0

        # read parameters
        with open(
            "parameters/parameters_test_combinedMomentEstimator.json", "r"
        ) as parameter_file:
            parameters = load(parameter_file)

        # build momentEstimator class
        test_me = CombinedMomentEstimator(**parameters["momentEstimatorInputDict"])

        # update power sums
        for i in range(len(list_values)):
            test_me.update([[list_values[i]]])

        # test
        xmc_estimation_order_1 = test_me.value(
            order=1, isCentral=True, isErrorEstimationRequested=False
        )
        xmc_estimation_order_2 = test_me.value(
            order=2, isCentral=True, isErrorEstimationRequested=False
        )
        xmc_error_order_1 = test_me.value(
            order=1, isCentral=True, isErrorEstimationRequested=True
        )
        self.assertEqual(xmc_estimation_order_1, true_estimation_order_1)
        self.assertEqual(xmc_estimation_order_2, true_estimation_order_2)
        self.assertEqual(xmc_error_order_1, true_error_order_1)


class TestMultiMomentEstimator(unittest.TestCase):
    """
    Class to test the class MultiMomentEstimator.

    Random testing is favoured: random numbers are drawn at every test execution, and sample
    estimations are compared with reference values deduced from knowledge of the probability
    distribution. Such test have a small probability of false negative (i.e. failure).

    For tests where this is not possible, samples are read from stored data and compared to
    reference values computed from these samples by means assumed correct.
    """

    @classmethod
    def setUpClass(cls):
        cls._randomGenerator = np.random.default_rng()
        cls.tolerance = 0.1

    def setUp(self):
        self.variableDimension = int(self._randomGenerator.uniform(1, 11, 1))
        self.mean = self._randomGenerator.uniform(-10, 10, self.variableDimension)
        self.variance = self._randomGenerator.uniform(0.1, 0.5, self.variableDimension)

        # Set number of samples
        # Accepted probability that estimated expectation does not satisfy tolerance
        failureProbability = 10 ** -4
        requiredStD = normalInverseCDF(1 - failureProbability / 2) / self.tolerance
        # Choose number of samples to get this failure probability
        self.numberOfSamples = int(np.max(self.variance * requiredStD ** 2))
        # Safety factor for higher-order moments
        self.numberOfSamples = self.numberOfSamples ** 2
        # Safety bounds
        self.numberOfSamples = max(10 ** 3, min(10 ** 6, self.numberOfSamples))

    def distance(self, x: List[Union[float, int]]) -> float:
        """
        Distance function used to quantify error.
        """
        return float(np.max(np.abs(x)))

    def test_estimation_random(self):
        """
        Randomised testing of the estimations. False failures are possible.
        """
        for dim, order, isError in itproduct((0, 1), (1, 2, 3, 4), (False, True)):
            not_implemented = dim == 1 or (isError and order > 2)
            if not_implemented:
                # Nothing to do
                continue
            if isError:
                reference = gaussianHStatVariance(self.variance, order, self.numberOfSamples)
            else:
                if order == 1:
                    # order 1: not actually a h-statistics
                    reference = gaussianRawMoment(self.mean, self.variance, order)
                else:
                    reference = gaussianCentralMoment(self.variance, order)
            me = MultiMomentEstimator(order=order)
            me.update(self._samples(dim))
            estimation = get_value_from_remote(me.multiValue(order, isError))
            # Test each component individually
            for c, (est, ref) in enumerate(zip(estimation, reference)):
                # Consider relative error if possible
                tol = abs(self.tolerance * ref)
                if tol == 0:
                    # Absolute error is considered
                    tol = self.tolerance
                with self.subTest(
                    msg=(
                        f"{'Variance of ' if isError else ''}{'Delta ' if dim==1 else ''}"
                        f"h-statistics of order {order}, component {c}"
                    ),
                    indexSetDimension=dim,
                    statisticalOrder=order,
                    errorEstimation=isError,
                    component=c,
                ):
                    self.assertAlmostEqual(est, ref, delta=tol)

    def test_estimation_deterministic(self):
        """
        Test estimation with respect to reference data.
        The reference data is generated by multi-moment_test_data.py
        """
        # Data for deterministic tests
        # The data is assumed to be small, so we store it all
        with open("parameters/multi-moment_test_data.json", "r") as f:
            referenceData = load(f)

        for dim, order, isError in itproduct((0, 1), (1, 2, 3, 4), (False, True)):
            referenceKey = f"{'Delta-' if dim == 1 else ''}h{order}{'_var' if isError else ''}"
            reference = referenceData[referenceKey]
            # Compute estimation
            estimator = MultiMomentEstimator(order=order)
            samples = referenceData["samples"]
            if dim == 0:
                # Extract samples from coarser (i.e. second) level, but preserve depth
                samples = [[s[1]] for s in samples]
            estimator.update(samples)
            estimation = get_value_from_remote(estimator.multiValue(order, isError))
            # Test each component individually
            for c, (est, ref) in enumerate(zip(estimation, reference)):
                if ref != 0:
                    # Consider relative error if possible
                    tol = abs(self.tolerance * ref)
                else:
                    # Absolute error is considered
                    tol = self.tolerance
                with self.subTest(
                    msg=(
                        f"{'Variance of ' if isError else ''}{'Delta ' if dim==1 else ''}"
                        f"h-statistics of order {order}, component {c}"
                    ),
                    indexSetDimension=dim,
                    statisticalOrder=order,
                    errorEstimation=isError,
                    component=c,
                ):
                    self.assertAlmostEqual(est, ref, delta=tol)

    def test_update(self):
        # Choose a reasonable number of events
        numberOfEvents = int(self._randomGenerator.uniform(1, 11, 1))
        # Test for an index set of dimension 0 (single-level) and 1 (multi-level)
        for dim in (0, 1):
            # Draw samples
            samples = self._samples(dim, numberOfEvents)
            # Compute reference power sums
            if dim == 0:
                referencePowerSums = powerSumsDimension0(samples)
            elif dim == 1:
                referencePowerSums = powerSumsDimension1(samples)
            referencePowerSums = get_value_from_remote(referencePowerSums)
            # Create test estimator and update it with samples
            estimator = MultiMomentEstimator(order=4)
            estimator.update(samples)
            # Get power sums to be tested
            testedPowerSums = get_value_from_remote(estimator._powerSums)
            # Run sub-test for the current index set dimension
            with self.subTest(
                msg=f"{'Multi' if dim==1 else 'Single'}-level update", indexSetDimension=dim
            ):
                # Check that the sample counter is right
                self.assertEqual(estimator.sampleNumber(), numberOfEvents)
                # Check every power sums against its reference value,
                # with a tolerance for tiny numerical errors
                for key, reference in referencePowerSums.items():
                    self.assertAlmostEqual(
                        self.distance(testedPowerSums[key] - reference), 0, delta=1 ** -15
                    )

    def test_construction_isParallel(self):
        """
        Tests the setting of attributes _isUpdateParallel and _isEstimationParallel of
        MultiMomentEstimator from keywords isParallel, isUpdateParallel and
        isEstimationParallel.
        """
        for is_parallel, is_update_parallel, is_estimation_parallel in itproduct(
            (None, False, True), repeat=3
        ):
            kwArgs = {
                "order": 1,
                "indexSetDimension": 1,
                "variableDimension": 2,
            }
            if is_parallel is not None:
                kwArgs["isParallel"] = is_parallel
            if is_update_parallel is not None:
                kwArgs["isUpdateParallel"] = is_update_parallel
            if is_estimation_parallel is not None:
                kwArgs["isEstimationParallel"] = is_estimation_parallel
            estimator = MultiMomentEstimator(**kwArgs)
            with self.subTest(
                isParallel=is_parallel,
                isUpdateParallel=is_update_parallel,
                msg=(
                    f"Testing construction with isParallel {is_parallel} "
                    f"and isUpdateParallel {is_update_parallel}"
                ),
            ):
                if is_update_parallel is not None:
                    self.assertEqual(estimator._isUpdateParallel, is_update_parallel)
                elif is_parallel is not None:
                    self.assertEqual(estimator._isUpdateParallel, is_parallel)
                else:
                    self.assertTrue(estimator._isUpdateParallel)
            with self.subTest(
                isParallel=is_parallel,
                isEstimationParallel=is_estimation_parallel,
                msg=(
                    f"Testing construction with isParallel {is_parallel} "
                    f"and isEstimationParallel {is_update_parallel}"
                ),
            ):
                if is_estimation_parallel is not None:
                    self.assertEqual(estimator._isEstimationParallel, is_estimation_parallel)
                elif is_parallel is not None:
                    self.assertEqual(estimator._isEstimationParallel, is_parallel)
                else:
                    self.assertTrue(estimator._isEstimationParallel)

    def pycompss_quake_test_isParallel(self):
        """
        Tests the parallel behaviours of estimation and update of MultiMomentEstimator,
        declared by its attributes _isUpdateParallel and _isEstimationParallel. This
        behaviour is either parallel or serial. Both options are tests successively in
        sub-tests. This test only makes sense when run in a parallel framework. Otherwise, the
        sub-test on parallel behaviour will always fail and the one on serial behaviour will
        always pass.

        """
        isd = 1  # index set dimension
        dims = (9, 2 ** isd, 2)  # events, levels, components
        commonArgs = {
            "order": 1,
            "indexSetDimension": isd,
            "variableDimension": dims[2],
        }
        # Test parallel behaviour
        for is_update_parallel, is_estimation_parallel in itproduct((False, True), repeat=2):
            # Samples have to be generated again every time or COMPSs fails to find them twice
            future_samples = futureArray(self._randomGenerator.normal(0, 1, dims).tolist())
            estimator = MultiMomentEstimator(
                isUpdateParallel=is_update_parallel,
                isEstimationParallel=is_estimation_parallel,
                **commonArgs,
            )
            with self.subTest(
                isUpdateParallel=is_update_parallel,
                msg=f"Testing {'parallel' if is_update_parallel else 'sequential'} update",
            ):
                estimator.update(future_samples)
                if is_update_parallel:
                    with self.assertRaises(
                        TypeError,
                        msg="MultiMomentEstimator._powerSums contains lists instead of Future.",
                    ):
                        estimator._powerSums["10"][0]
                else:
                    self.assertIsInstance(estimator._powerSums["10"][0], float)
            with self.subTest(
                isEstimationParallel=is_estimation_parallel,
                isUpdateParallel=is_update_parallel,
                msg=(
                    f"Testing {'parallel' if is_update_parallel else 'sequential'} "
                    "estimation, with "
                    f"{'parallel' if is_update_parallel else 'sequential'} update"
                ),
            ):
                v = estimator.value(1)
                if is_estimation_parallel:
                    with self.assertRaises(
                        TypeError,
                        msg=(
                            "MultiMomentEstimator.value returned "
                            f"{type(v)} instead of Future."
                        ),
                    ):
                        v += 1
                else:
                    self.assertIsInstance(v, float)

    def _samples(self, dimension: int, number: Union[int, None] = None) -> SampleArray:
        """
        Generates the requested number of independant samples,
        using the random generator of the object.
        """
        if number is None:
            number = self.numberOfSamples
        dims = (number, 2 ** dimension, self.variableDimension)
        # TODO Fix: this gives uncorrelated samples for an index set of dimension 1
        samples = self._randomGenerator.normal(self.mean, np.sqrt(self.variance), dims)
        # Ensure output type
        return samples.tolist()


class TestMultiCombinedMomentEstimator(unittest.TestCase):
    """
    Test case for the class MultiCombinedMomentEstimator.
    """

    @classmethod
    def setUpClass(cls):
        cls.order = 4
        cls.pointPerEvent = 10
        cls.eventNumber = 6
        # Shape of sample array: (random events, MC levels, variable components, power sums, 1)
        # This is as returned by KratosSolverWrapper.solve (except for counter, see below)
        power_sums_per_event = np.tile(
            [[[20.0], [40.0]], [[25.0], [45.0]]], (cls.eventNumber, 1, 1, cls.order, 1)
        )
        cls.variableDimension = power_sums_per_event.shape[2]
        # Add counter of number of terms in power sums for each component,
        # and store in class attribute
        cls.samples = power_sums_per_event.tolist()
        for perEvent in cls.samples:
            for perLevel in perEvent:
                for perComponent in perLevel:
                    perComponent.append(cls.pointPerEvent)
        # Reference values (deterministic)
        # First we remove axes of length 1: 'MC levels' (2nd), and the last one (pointless?)
        reference_power_sums = power_sums_per_event.squeeze(axis=(1, -1))
        # We sum power sums across events
        reference_power_sums = np.sum(reference_power_sums, axis=0)
        # We transpose axes from (component, power) to (power, component)
        reference_power_sums = reference_power_sums.transpose()
        cls.referencePowerSums = {}
        for i, v in enumerate(reference_power_sums):
            cls.referencePowerSums[str(i + 1)] = v
        cls.referenceStatistics = {
            "h1": [2.0, 2.5],
            "h1_var": [0.0, -0.02966102],
            "h2": [0.0, -1.77966102],
            "h2_var": [0.64283737, 0.49816228],
            "h3": [-6.31209819, 0.0],
            "h3_var": [624.86742843, 1629.54759381],
            "h4": [38.57024209, 32.94956781],
            "h4_var": [-3624.77860536, 24223.33392998],
        }
        # TODO: these statistics are obviously absurd, because the power sums are

    def setUp(self):
        # Build estimator
        self._estimator = MultiCombinedMomentEstimator(
            order=self.order, variableDimension=self.variableDimension
        )
        # update power sums
        self._estimator.update(deepcopy(self.samples))

    def test_updateD0(self):
        # test sample number
        sample_count = get_value_from_remote(self._estimator.sampleNumber())
        self.assertEqual(sample_count, self.pointPerEvent * self.eventNumber)

        # test power sums
        power_sums = get_value_from_remote(self._estimator._powerSums)
        for key, ref in self.referencePowerSums.items():
            for c in range(self.variableDimension):
                with self.subTest(
                    msg=f"Single-level update, for component {c}", powerSum=key, component=c
                ):
                    self.assertEqual(power_sums[key][c], ref[c])

    def test_estimationD0(self):
        for order, isError in itproduct(range(1, self.order + 1), (False, True)):
            key = f"h{order}{'_var' if isError else ''}"
            estimation = get_value_from_remote(self._estimator.multiValue(order, isError))
            for c in range(self.variableDimension):
                with self.subTest(
                    msg=(
                        f"{'Variance of ' if isError else ''}h-statistics of order {order}, "
                        f"component {c}"
                    ),
                    powerSum=key,
                    component=c,
                ):
                    self.assertAlmostEqual(
                        estimation[c],
                        self.referenceStatistics[key][c],
                    )


def gaussianCentralMoment(
    variance: Union[float, List[float]], order: int
) -> Union[float, List[float]]:
    """
    Returns the central moment of any strictly positive order of a random variable
    following a normal distribution of given variance. Its mean matters not.
    Reference: https://math.stackexchange.com/a/92650
    """
    if order % 2 == 1:
        # Symmetric distribution: moments of odd orders are zero
        return 0 * variance
    else:
        return doubleFactorial(order - 1) * variance ** (order / 2)


def gaussianHStatVariance(
    gaussianVariance: float, hStatOrder: int, numberOfSamples: int
) -> float:
    """
    Returns the variance of the h-statistics for the given number of samples.
    The random variable is assumed to follow a normal distribution.
    """
    if hStatOrder == 1:
        # Not a h-statistics, but it is convenient to support it here
        return gaussianVariance / numberOfSamples
    elif hStatOrder == 2:
        gaussianM4 = gaussianCentralMoment(gaussianVariance, 4)
        return gaussianM4 / numberOfSamples - (gaussianVariance ** 2) * (
            numberOfSamples - 3
        ) / (numberOfSamples * (numberOfSamples - 1))
    else:
        raise ValueError(f"Not implemented for order {hStatOrder} (>2)")


def gaussianRawMoment(mean: float, variance: float, order: int) -> float:
    """
    Returns the raw moment of any strictly positive order of a random variable
    following a normal distribution of given mean and variance.
    Reference: https://math.stackexchange.com/a/1945494
    """
    moment = 0
    for i in range(0, order + 1):
        if order - i % 2 == 0:
            t1 = (order - i) // 2
            t2 = factorial(2 * t1) / (factorial(t1) * 2 ** t1)
            moment += binomial(order, i) * mean ** i * variance ** (t1) * t2
    return moment


@task(keep=True, returns=1, values={Type: COLLECTION_IN, Depth: 3})
def powerSumsDimension0(values: SampleArray) -> PowerSumsDict:
    """
    Produces reference values of power sums from samples corresponding to a Monte Carlo index
    of dimension 0.
    """
    v = np.array(values)
    # Power sums (initially an empty dictionary)
    ps = {}
    # Computations of power sums
    # Ugly, but this is the reference
    ps["1"] = np.sum(v, axis=0)
    ps["2"] = np.sum(v ** 2, axis=0)
    ps["3"] = np.sum(v ** 3, axis=0)
    ps["4"] = np.sum(v ** 4, axis=0)
    ps["5"] = np.sum(v ** 5, axis=0)
    ps["6"] = np.sum(v ** 6, axis=0)
    ps["7"] = np.sum(v ** 7, axis=0)
    ps["8"] = np.sum(v ** 8, axis=0)
    return ps


@task(keep=True, returns=1, values={Type: COLLECTION_IN, Depth: 3})
def powerSumsDimension1(values: SampleArray) -> PowerSumsDict:
    """
    Produces reference values of power sums from samples corresponding to a Monte Carlo index
    of dimension 1.
    """
    values = np.array(values)
    # Sum of values
    vs = values[:, 0] + values[:, 1]
    # Difference of values
    vd = values[:, 0] - values[:, 1]
    # Power sums (initially an empty dictionary)
    ps = {}
    # Computations of power sums
    # Ugly, but this is the reference
    ps["10"] = np.sum(vs, axis=0)
    ps["20"] = np.sum(vs ** 2, axis=0)
    ps["30"] = np.sum(vs ** 3, axis=0)
    ps["40"] = np.sum(vs ** 4, axis=0)
    ps["50"] = np.sum(vs ** 5, axis=0)
    ps["60"] = np.sum(vs ** 6, axis=0)
    ps["70"] = np.sum(vs ** 7, axis=0)
    ps["80"] = np.sum(vs ** 8, axis=0)
    ps["01"] = np.sum(vd, axis=0)
    ps["11"] = np.sum(vs * vd, axis=0)
    ps["21"] = np.sum(vs ** 2 * vd, axis=0)
    ps["31"] = np.sum(vs ** 3 * vd, axis=0)
    ps["41"] = np.sum(vs ** 4 * vd, axis=0)
    ps["51"] = np.sum(vs ** 5 * vd, axis=0)
    ps["61"] = np.sum(vs ** 6 * vd, axis=0)
    ps["71"] = np.sum(vs ** 7 * vd, axis=0)
    ps["02"] = np.sum(vd ** 2, axis=0)
    ps["12"] = np.sum(vs * vd ** 2, axis=0)
    ps["22"] = np.sum(vs ** 2 * vd ** 2, axis=0)
    ps["32"] = np.sum(vs ** 3 * vd ** 2, axis=0)
    ps["42"] = np.sum(vs ** 4 * vd ** 2, axis=0)
    ps["52"] = np.sum(vs ** 5 * vd ** 2, axis=0)
    ps["62"] = np.sum(vs ** 6 * vd ** 2, axis=0)
    ps["03"] = np.sum(vd ** 3, axis=0)
    ps["13"] = np.sum(vs * vd ** 3, axis=0)
    ps["23"] = np.sum(vs ** 2 * vd ** 3, axis=0)
    ps["33"] = np.sum(vs ** 3 * vd ** 3, axis=0)
    ps["43"] = np.sum(vs ** 4 * vd ** 3, axis=0)
    ps["53"] = np.sum(vs ** 5 * vd ** 3, axis=0)
    ps["04"] = np.sum(vd ** 4, axis=0)
    ps["14"] = np.sum(vs * vd ** 4, axis=0)
    ps["24"] = np.sum(vs ** 2 * vd ** 4, axis=0)
    ps["34"] = np.sum(vs ** 3 * vd ** 4, axis=0)
    ps["44"] = np.sum(vs ** 4 * vd ** 4, axis=0)
    ps["05"] = np.sum(vd ** 5, axis=0)
    ps["15"] = np.sum(vs * vd ** 5, axis=0)
    ps["25"] = np.sum(vs ** 2 * vd ** 5, axis=0)
    ps["35"] = np.sum(vs ** 3 * vd ** 5, axis=0)
    ps["06"] = np.sum(vd ** 6, axis=0)
    ps["16"] = np.sum(vs * vd ** 6, axis=0)
    ps["26"] = np.sum(vs ** 2 * vd ** 6, axis=0)
    ps["07"] = np.sum(vd ** 7, axis=0)
    ps["17"] = np.sum(vs * vd ** 7, axis=0)
    ps["08"] = np.sum(vd ** 8, axis=0)
    return ps


def futureArray(array: SampleArray):
    for i, j, k in itproduct(range(len(array)), range(len(array[0])), range(len(array[0][0]))):
        array[i][j][k] = returnInput_Task(array[i][j][k])
    return array
