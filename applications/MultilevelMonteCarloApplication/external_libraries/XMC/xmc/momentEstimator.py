from itertools import product as itproduct
from numpy import ceil, log2, zeros, float64

# Import xmc classes
from xmc.statisticalEstimator import StatisticalEstimator
from xmc.tools import dynamicImport
from exaqute import get_value_from_remote
from xmc.methodDefs_momentEstimator.types import *
from xmc.methodDefs_momentEstimator import (
    hStatistics,
    powerSums,
    updatePowerSums,
    updateCombinedPowerSums,
    computeCentralMoments,
    computeErrorEstimation,
)


class MomentEstimator(StatisticalEstimator):
    """
    This class estimates raw and central moments up to a given order (including
    expectations, obviously). These estimations are computed using h-statistics,
    so we store (and update) the power sums from which they are computed.
    It inherits from StatisticalEstimator; the description below
    describes only its differences with it.

    Description of attributes
    _________________________

    * order: integer.
        Mandatory.
        All statistical moments up to this order are estimated. It determines the power of the sums of realisations to be computed; see MomentEstimator.powerSumsOrder.

    * indexSetDimension: integer.
        Local dimension of the Monte Carlo index set, i.e. the number of solvers sending samples to this estimator is 2^self.indexSetDimenion. It is only used to initialise self.powerSums and define the default value of self._updatedPowerSums.

    * powerSums: List[float]
        List of power sums required to compute the moments.

    * _updatedPowerSums : callable.
        Default: xmc.methodDefs_momentEstimator.updatePowerSums.updatePowerSumsOrder{o}Dimension{d}, where o is the order and d is the indexSetDimension. It may be set by the user with the key \"updatedPowerSums\".
        Method which updates power sums with additional sample values; called by self.update (see also its documentation).

    * _centralMomentComputer : callable.
        Default: xmc.methodDefs_momentEstimator.computeCentralMoments.centralMomentWrapper. It may be set by the user with the key \"centralMomentComputer\".
        Method which computes the desired statistics (central moments) from the available power sums; called by self.value (see also its documentation).

    * _centralMomentErrorComputer : callable.
        Default: xmc.methodDefs_momentEstimator.computeErrorEstimation.centralMomentErrorWrapper. It may be set by the user with the key \"centralMomentErrorComputer\".
        Method which computes the error of the desired statistics (central moments) from the available power sums; called by self.value (see also its documentation).

    * _rawMomentComputer : callable.
        Default: None. It may be set by the user with the key \"rawMomentComputer\".
        Method which computes the desired statistics (raw moments) from the available power sums; called by self.value (see also its documentation).

    * _rawMomentErrorComputer : callable.
        Default: None. It may be set by the user with the key \"rawMomentErrorComputer\".
        Method which computes the error of the desired statistics (raw moments) from the available power sums; called by self.value (see also its documentation).
    """

    def __init__(self, **keywordArgs):
        """
        Keyword arguments required:
            * order
            * indexSetDimension

        Optional keyword arguments:
            * updatedPowerSums
            * centralMomentComputer
            * centralMomentErrorComputer
            * rawMomentComputer
            * rawMomentErrorComputer
        """
        # Parent constructor
        super().__init__(**keywordArgs)

        # Attributes
        # maximum order to which the moments are computed
        self.order = keywordArgs.get("order")
        # TODO indexSetDimension is only used to initialise self.powerSums
        # and define the default value of self._updatedPowerSums.
        # Should it be a required input argument?
        indexSetDimension = keywordArgs.get("indexSetDimension")
        self.indexSetDimension = indexSetDimension
        # list of power sums required to compute the moments
        # TODO - This is a temporary fix until when sample_moments.py is interfaced
        if indexSetDimension == 0:
            self.powerSums = [[None] for _ in range(self.powerSumsOrder())]
        elif indexSetDimension == 1:
            self.powerSums = [
                [None for _ in range(i + 2)] for i in range(self.powerSumsOrder())
            ]
        else:
            self.powerSums = None

        # Methods
        basePath = "xmc.methodDefs_momentEstimator."
        self._centralMomentComputer = dynamicImport(
            keywordArgs.get(
                "centralMomentComputer",
                basePath + "computeCentralMoments.centralMomentWrapper",
            )
        )
        self._centralMomentErrorComputer = dynamicImport(
            keywordArgs.get(
                "centralMomentErrorComputer",
                basePath + "computeErrorEstimation.centralMomentErrorWrapper",
            )
        )
        self._rawMomentComputer = dynamicImport(keywordArgs.get("rawMomentComputer", None))
        self._rawMomentErrorComputer = dynamicImport(
            keywordArgs.get("rawMomentErrorComputer", None)
        )
        # updatedPowerSums method (conditional default value)
        defaultUpdater = (
            basePath + "updatePowerSums.updatePowerSums" "Order{o}Dimension{d}"
        ).format(o=self.powerSumsOrder(), d=self.indexSetDimension)
        self._updatedPowerSums = dynamicImport(
            keywordArgs.get("updatedPowerSums", defaultUpdater)
        )

    def update(self, samples: SampleArray):
        """
        Method updating the power sums given new samples.

        Inputs:
        - self: an instance of the class
        - samples: list containing new samples (dimension at least 3)
        """

        dimension = log2(len(samples[0]))  # len(power_sum[0])=2**dimension
        # Let's unpack the nested list self.powerSums
        # First all elements in self.powerSums[0], then those in self.powerSums[1], etc.
        # Idem for samples
        # TODO is it the order in which we want these lists?
        power_sums = [item for sublist in self.powerSums for item in sublist]
        power_sums = self._updatedPowerSums(samples, *power_sums)
        power_sums = list(power_sums)
        for i in range(len(self.powerSums)):
            for j in range(len(self.powerSums[i])):
                self.powerSums[i][j] = power_sums.pop(0)
        self._sampleCounter = self._sampleCounter + len(samples)

    def value(self, order, isCentral=None, isErrorEstimationRequested=False):
        """
        Method returning a specific order moment, its error estimate or both.

        Inputs:
        - self: an instance of the class
        - order: order of the moment
        - isCentral: setting if the working moment is central or raw
        - computeErrorEstimation: setting if computing or not the errorEstimation
        """

        if isCentral is None:
            isCentral = order > 1
        # TODO the [] around the returns in the following is to ensure consistency with
        # the inputs expected by the rest of the program. This is becuase the quantity
        # returned can be either a scalar estimation value or a list of values of function
        # evalutaions at certain points in a domain (eg. CDF)
        if isCentral:
            return self._computeCentralMoment(order, isErrorEstimationRequested)
        else:
            return self._computeRawMoment(order, isErrorEstimationRequested)

    def errorEstimation(self, order, isCentral):
        """
        Method returning the error estimation of working moment from stored power sums.

        Inputs:
        - self: an instance of the class
        - order: order of the moment
        """

        # TODO not normally used any more. Keep ?
        # TODO - before inferring dimension from powerSums, need to ensure that
        # update method has been run at least a certain number of times required
        # to compute the moment of requested order
        dimension = log2(len(self.powerSums[0]))  # len(power_sum[0])=2**dimension
        # Let's unpack the nested list self.powerSums
        # First all elements in self.powerSums[0], then those in self.powerSums[1], etc.
        # TODO is it the order in which we want this list?
        power_sums = [item for sublist in self.powerSums for item in sublist]
        if isCentral:
            return self._centralMomentErrorComputer(
                dimension, order, *[*power_sums, self.sampleNumber()]
            )
        else:
            return self._rawMomentErrorComputer(
                dimension, order, *[*power_sums, self.sampleNumber()]
            )

    def _computeCentralMoment(self, order, isErrorEstimationRequested):
        """
        Method returning the central moment of working moment from stored power sums.

        Inputs:
        - self: an instance of the class
        - order: order of the moment
        """

        dimension = log2(len(self.powerSums[0]))  # len(power_sum[0])=2**dimension
        # Let's unpack the nested list self.powerSums
        # First all elements in self.powerSums[0], then those in self.powerSums[1], etc.
        # TODO is it the order in which we want this list?
        if isErrorEstimationRequested:
            power_sums = [item for sublist in self.powerSums[: 2 * order] for item in sublist]
            output = self._centralMomentErrorComputer(
                dimension, order, *[*power_sums, self.sampleNumber()]
            )
        else:
            power_sums = [item for sublist in self.powerSums[:order] for item in sublist]
            output = self._centralMomentComputer(
                dimension, order, *[*power_sums, self.sampleNumber()]
            )
        return output

    def _computeRawMoment(self, order, isErrorEstimationRequested):
        """
        Compute raw moment of requested order from stored power sums.
        """

        dimension = log2(len(self.powerSums[0]))  # len(power_sum[0])=2**dimension
        # Let's unpack the nested list self.powerSums
        # First all elements in self.powerSums[0], then those in self.powerSums[1], etc.
        # TODO is it the order in which we want this list?
        power_sums = [item for sublist in self.powerSums for item in sublist]
        output = self._rawMomentComputer(dimension, order, *power_sums, self.sampleNumber())
        if isErrorEstimationRequested:
            output = [
                output,
                self._rawMomentErrorComputer(
                    dimension, order, *power_sums, self.sampleNumber()
                ),
            ]
        return output

    def powerSumsOrder(self):
        """
        Returns the maximum order to which the power sums are stored to compute estimation and
        error of moments up to order self.order.
        """

        return 2 * self.order

    def reset(self):
        """
        Reset the value of power sums to zero. Reset the number of samples
        """

        noneList = [[None] * len(powerSum) for powerSum in self.powerSums]
        self.powerSums = noneList
        self._sampleCounter = 0


class CombinedMomentEstimator(MomentEstimator):
    """
    This class estimates raw and central moments up to a given order (including expectations, obviously).
    These estimations are computed using h-statistics, so we store (and update) the power sums from which they are computed.
    The computed statistics are the composition of sampling and time statistical processes.
    The first difference with respect to MomentEstimator is the update method, which is supposed to receive power sums instead of sample values.
    Another difference with respect to MomentEstimator is that, for multi-index moment estimators (as Multilevel Monte Carlo),
    we can compute only power sums of order a of the form (Sa0,S0a), which requires to estimate h-statistics
    exploiting equations as eq. (4) of [S. Krumscheid et al. / Journal ofComputational Physics 414 (2020) 109466].

    It inherits from MomentEstimator; the description below
    describes only its differences with it.

    * _updatedPowerSums : callable.
        Default: xmc.methodDefs_momentEstimator.updateCombinedPowerSums.updatePowerSumsOrder{o}Dimension{d}, where o is the order and d is the indexSetDimension. It may be set by the user with the key \"updatedPowerSums\".
        Method which updates power sums with additional sample values; called by self.update (see also its documentation).

    * _centralMomentComputer : callable.
        Default: xmc.methodDefs_momentEstimator.computeCentralMoments.centralCombinedMomentErrorWrapper. It may be set by the user with the key \"centralMomentComputer\".
        Method which computes the desired statistics (central moments) from the available power sums; called by self.value (see also its documentation).

    * _centralMomentErrorComputer : callable.
        Default: xmc.methodDefs_momentEstimator.computeErrorEstimation.centralMomentErrorWrapper. It may be set by the user with the key \"centralMomentErrorComputer\".
        Method which computes the error of the desired statistics (central moments) from the available power sums; called by self.value (see also its documentation).

    * _rawMomentComputer : callable.
        Default: None. It may be set by the user with the key \"rawMomentComputer\".
        Method which computes the desired statistics (raw moments) from the available power sums; called by self.value (see also its documentation).

    * _rawMomentErrorComputer : callable.
        Default: None. It may be set by the user with the key \"rawMomentErrorComputer\".
        Method which computes the error of the desired statistics (raw moments) from the available power sums; called by self.value (see also its documentation).

    Methods:
    - update: method updating power sums and number of realizations.
    """

    def __init__(self, **keywordArgs):
        # Parent constructor
        super(CombinedMomentEstimator, self).__init__(**keywordArgs)

        if self.indexSetDimension == 0:
            self.powerSums = [[None] for _ in range(self.powerSumsOrder())]
        elif self.indexSetDimension == 1:
            self.powerSums = [
                [None for _ in range(0, 2)] for _ in range(self.powerSumsOrder())
            ]
        else:
            self.powerSums = None

        # Methods
        basePath = "xmc.methodDefs_momentEstimator."
        self._centralMomentComputer = dynamicImport(
            keywordArgs.get(
                "centralMomentComputer",
                basePath + "computeCentralMoments.centralCombinedMomentWrapper",
            )
        )
        self._centralMomentErrorComputer = dynamicImport(
            keywordArgs.get(
                "centralMomentErrorComputer",
                basePath + "computeErrorEstimation.centralCombinedMomentErrorWrapper",
            )
        )
        self._rawMomentComputer = dynamicImport(keywordArgs.get("rawMomentComputer", None))
        self._rawMomentErrorComputer = dynamicImport(
            keywordArgs.get("rawMomentErrorComputer", None)
        )
        # updatedPowerSums method (conditional default value)
        defaultUpdater = (
            basePath + "updateCombinedPowerSums.updatePowerSums" "Order{o}Dimension{d}"
        ).format(o=self.powerSumsOrder(), d=self.indexSetDimension)
        self._updatedPowerSums = dynamicImport(
            keywordArgs.get("updatedPowerSums", defaultUpdater)
        )

    def update(self, samples):
        """
        Method updating the power sums and number of realizations, given new samples. The new
        data is passed through as a list containing power sums, up to a given order, and the
        number of samples from which it was computed. For example, let us discretise a scalar
        time signal u into M time steps to obtain a vector U of length M. A power sum of order
        a is computed as S_a = \sum_{i=1}^{M} U[i]^a.
        Samples will have the following shape:
        [[[[S1], [S2], M]]] ,
        where S1 and S2 are power sums of order one and two, respectively, and M is the total
        number of terms in the power sum (equivalent to number of time steps times realizations).

        Inputs:
        - self: an instance of the class
        - samples: list containing new sample power sums.
        """

        # Let's unpack the nested list self.powerSums
        # First all elements in self.powerSums[0], then those in self.powerSums[1], etc.
        # Idem for samples
        power_sums = [item for sublist in self.powerSums for item in sublist]
        sample_counter_power_sums = self._updatedPowerSums(
            self._sampleCounter, samples, *power_sums
        )
        sample_counter = sample_counter_power_sums[0]
        power_sums = list(sample_counter_power_sums[1:])
        for i in range(len(self.powerSums)):
            for j in range(len(self.powerSums[i])):
                self.powerSums[i][j] = power_sums.pop(0)
        self._sampleCounter = sample_counter


class MultiMomentEstimator(StatisticalEstimator):
    """
    This class estimates statistical moments of multi-valued real random variables, from
    independent realisations. It inherits from StatisticalEstimator; the description below
    describes only its differences with it.

    Description of attributes
    _________________________

    * order : int
    Mandatory.
    All statistical moments up to this order are estimated. It determines the power of the sums
    of realisations to be computed; see MultiMomentEstimator.powerSumsOrder.

    * _referenceComponent : int
    Default: 0.
    MultiMomentEstimator.value must return a float. Of all the components of the multi-valued
    random variable, MultiMomentEstimator.value will return the requested statistic of the one
    of index MultiMomentEstimator._referenceComponent.

    * _isEstimationParallel : bool
    Default: True.
    Whether the estimation process is to run in parallel. If not, power sums are synchronised
    synchronised before computing the estimation.

    * _isEstimationParallel : bool
    Default: True.
    Whether the update process is to run in parallel. If not, samples are synchronised before
    the update.

    * _variableDimension : int
    Mandatory if self._isParallel; otherwise, automatically assigned.  The dimension of the
    image space (i.e. number of scalar components) of the multi-valued random variable. If
    self._isParallel, then it must be supplied to the constructor; if not, it can be guessed
    from the first samples.

    * _indexSetDimension : int
    Automatically assigned.
    Local dimension of the Monte Carlo index set, i.e. the number of solvers sending samples to
    this estimator is 2^self._indexSetDimenion. It is not supplied at construction but deduced
    from the first samples.

    * _powerSums : DefaultDict[str, List[float]]
    Default: see self._initialisePowerSums.
    This dictionary stores sums of powers of sample values, and can be supplied at
    construction. The summation and exponentiation are component-wise, so that the values in
    each entry is an array of dimension self._variableDimension. The keys are strings of
    integers, of length 2^self._indexSetDimension; each integer is a power of a multivariate
    sum or difference of multi-valued random variables across solver levels.
    Examples:
    - if self._indexSetDimension is 0: there is one solver level (multi-valued random variable
    X), so the keys are of the form 'p' where p is the integer power of the sums stored in this
    entry.
    self._powerSums['p'] = sum_{e in events} X(e)^p
    - if self._indexSetDimension is 1: there are two solver levels (multi-valued random
    variables X and Y), so the keys are of the form 'pq' where p and q are respectively the
    integer powers of the sum and difference of multi-valued random variables across solver
    levels.
    self._powerSums['pq'] = sum_{e in events} (X(e)+Y(e))^p (X(e)-Y(e))^q

    * _powerSumsUpdater : callable
    Default: xmc.methodDefs_momentEstimator.powerSums.updatedPowerSums.
    Function which updates power sums with additional sample values; called by self.update (see
    also its documentation). It must abide by the following specifications; see default
    function for an example.
    Input arguments:
    1. dictionary of power sums: DefaultDict[str, List[float]]
    2. array of samples: List[List[List[float]]]
    Output arguments:
    1. dictionary of power sums (updated): DefaultDict[str, List[float]]

    * _estimation: callable
    Default: xmc.methodDefs_momentEstimator.hStatistics.hStatistics.
    Function which computes the desired statistics from the available power sums; called by
    self.multiValue (see also its documentation). It must abide by the following
    specifications; see default function for an example.
    Input arguments:
    1. components of the multi-valued random variable to consider: Union[int, slice,
    Iterable[int]]
    2. dictionary of power sums: DefaultDict[str, List[float]].
    3. number of random events having contributed to the power sums: int.
    4. statistical order of the moment estimated: int.
    5. whether the error on the estimation is requested, instead of the estimation: bool.
    6. whether the moment estimated is central or raw: bool.
    Output arguments:
    1. statistics for every component considered: Union[List[float],float] (float if single
    component)

    """

    def __init__(self, **keywordArgs):
        """
        Keyword arguments required:
        - order

        Keyword arguments required for parallel update:
        - variableDimension

        Optional keyword arguments:
        - referenceComponent
        - isUpdateParallel
        - isEstimationParallel
        - isParallel
        - estimationFunction
        - powerSumsUpdater
        - powerSums

        Notes:
        - isParallel: value is assigned to both self._isUpdateParallel and self._isEstimationParallel. If the keyword arguments isUpdateParallel or isEstimationParallel are given, they take priority.
        """
        # Parent constructor
        super().__init__(**keywordArgs)
        #
        # Attributes
        # maximum order to which the moments are computed
        self.order = keywordArgs.get("order")
        # Index of reference component for the algorithm: error estimation and so on
        self._referenceComponent = keywordArgs.get("referenceComponent", 0)
        # Whether to run computations in parallel
        isParallel = keywordArgs.get("isParallel", True)
        self._isUpdateParallel = keywordArgs.get("isUpdateParallel", isParallel)
        self._isEstimationParallel = keywordArgs.get("isEstimationParallel", isParallel)
        # The following two attributes will be initialised upon first update,
        # but if the samples are Future objects, the variable dimension must be given
        self._variableDimension = keywordArgs.get("variableDimension", None)
        self._indexSetDimension = None
        #
        # Methods
        powerSumsUpdaterName = keywordArgs.get(
            "powerSumsUpdater", "xmc.methodDefs_momentEstimator.powerSums.updatedPowerSums"
        )
        # Modules for methods to compute moments
        estimationName = keywordArgs.get(
            "estimationFunction", "xmc.methodDefs_momentEstimator.hStatistics.hStatistics"
        )
        #
        # Non-empty initialisation
        # Enforce parallelisation
        if self._isUpdateParallel and powerSumsUpdaterName[-6:-1] != "_Task":
            powerSumsUpdaterName += "_Task"
        if self._isEstimationParallel and estimationName[-6:-1] != "_Task":
            estimationName += "_Task"
        self._powerSumsUpdater = dynamicImport(powerSumsUpdaterName)
        self._estimation = dynamicImport(estimationName)
        # list of power sums used to compute the moments
        self._powerSums = keywordArgs.get("powerSums", None)
        # in case power sums are input
        if self._powerSums is not None:
            self._indexSetDimension = self._guessIndexSetDimension()
            self._variableDimension = self._guessVariableDimension()

    def _initialisePowerSums(self):
        # Generate all combinations of self._indexSetDimension+1 integers between 0
        # and twice the maximal order of moments to compute.
        maxPower = 2 * self.order
        powers = itproduct(range(maxPower + 1), repeat=2 ** self._indexSetDimension)
        powerSumsKeys = [
            # concatenate as strings each list of integers
            "".join(map(str, p))
            # from every combination of integers whose sum is in ]0,maxPower]
            for p in powers
            if 0 < sum(p) <= maxPower
        ]
        # Generate list of dictionaries from these keys;
        # Initial value: NumPy array of zeros, in long format to avoid overflow.
        self._powerSums = {
            k: zeros(self._variableDimension, dtype=float64) for k in powerSumsKeys
        }

    def _guessIndexSetDimension(self, singleSample: List[List[float]] = None) -> int:
        """
        Return the dimension of the index set considered by this estimator. If unknown, returns
        None. Relies on the keys of dictionary self._powerSums being of length dimension+1

        Alternatively: if input argument is provided, guess index set dimension from it.

        Input argument:
        - singleSample (optional): set of correlated samples (list) for a single event.
        Expected to be of length 2^dimension

        Output argument: integer, or None if no dimension could be computed.
        """
        if singleSample is not None:
            # This list should be of length 2^dimension
            return int(ceil(log2(len(singleSample))))
        if self._powerSums is not None:
            # Use the length of first key of dictionary self._powerSums
            # It should be of length 2^dimension
            return int(ceil(log2(len(next(iter(self._powerSums))))))
        # If we reach here, then there is no information to process
        return None

    def _guessVariableDimension(self, singleSample: List[List[float]] = None) -> int:
        """
        Return the dimension of the space of values of the multi-valued random variable.

        Input argument:
        - singleSample (optional): set of correlated samples (list) for a single event.
        Every element of that list is expected to be a list of length equal to the dimension of
        the multi-valued random variable.

        Output argument: integer, or None if no dimension could be computed

        An error will be raised if the sample is a Future object. In that case, the variable
        dimension must be given. Use the variableDimension keyword at instantiation..
        """
        if singleSample is not None:
            # This list should contain sublists of length equal to the variable dimension
            if not isinstance(singleSample, list) or not isinstance(singleSample[0], list):
                raise TypeError(f"Expected a list of lists; received a {type(singleSample)}.")
            return len(singleSample[0])
        if self._powerSums is not None:
            # Get any key (first, here)
            anyKey = next(iter(self._powerSums))
            # The value of this key must be a list. Return its length.
            return len(self._powerSums[anyKey])
        # If we reach here, then there is no information to process
        return None

    def update(self, samples: SampleArray):
        """
        Method to update the power sums from new samples.

        Input arguments:
        1. samples: nested list of depth ≥ 3. (S_ijk) with: i the random event; j the solver
        level; k the component of the multi-valued random variable.
        """
        # Synchronise if this estimator does support parallel computations
        if not self._isUpdateParallel:
            samples = get_value_from_remote(samples)
        # Check proper format
        if not isinstance(samples, list) or not isinstance(samples[0], list):
            raise TypeError("Input argument is expected to be a list of lists")
        # At first update, initialise if necessary
        if self.sampleNumber() < 1:
            # Guess missing information from sample set
            if self._indexSetDimension is None:
                self._indexSetDimension = self._guessIndexSetDimension(samples[0])
            if self._variableDimension is None:
                self._variableDimension = self._guessVariableDimension(samples[0])
            # Initialise power sums and their updater if they are not set
            if self._powerSums is None:
                self._initialisePowerSums()
        # Proceed to update
        self._powerSums = self._powerSumsUpdater(self._powerSums, samples)
        self._addData(samples)

    def value(
        self, order: int, isError: bool = False, isCentral: Union[bool, None] = None
    ) -> float:
        """
        Function returning the current estimation of a statistical moment, or its error estimation.
        This function only returns estimations for the component indicated by self._referenceComponent.

        Inputs: same as MultiMomentEstimator.multiValue, except for the 'component',
        which is not an input argument of MultiMomentEstimator.value.
        """
        return self.multiValue(order, isError, self._referenceComponent, isCentral)

    def multiValue(
        self,
        order: int,
        isError: bool = False,
        component: ListIndex = slice(None, None),
        isCentral: Union[bool, None] = None,
    ) -> HStatistics:
        """
        Function returning a specific order moment, its error estimate or both.

        Inputs:
        1. order: order of the moment
        2. isError: whether to return the error of the estimation (as opposed to the estimation
        itself).
        3. component: index (or indices) of components of the multi-valued random variable to
        consider
        4. isCentral: whether the moment is central (as opposed to raw)

        Outputs: float if component is an integer, list of floats otherwise.
        """
        # Assign missing values
        if isCentral is None:
            # Assume that any moment of order > 1 is requested as central
            isCentral = order > 1
        #
        if not self._isEstimationParallel:
            self._powerSums = get_value_from_remote(self._powerSums)
        return self._estimation(
            component,
            self._powerSums,
            self.sampleNumber(),
            order,
            isError,
            isCentral,
        )


class MultiCombinedMomentEstimator(MultiMomentEstimator):
    """
    This class estimates contributions to statistical 'combined' moments of multi-valued real
    random variables, from independent realisations. These contribution take the form of
    h-statistics, or differences of such. They are computed from monovariate power sums. The
    computed statistics are the composition of sampling and time statistical processes. This
    class inherits from MultiMomentEstimator; the description below describes only its
    differences with it.

    Description of attributes
    _________________________

    * _powerSums: Union[DefaultDict[str, PowerSumsDict], PowerSumsDict] with
    PowerSumsDict = DefaultDict[str, Union[float, List[float]]]
    See the documentation of this attribute of the parent class. The only difference is if
    self._indexSetDimension is 1: since only monovariate power sums are used, self._powerSums
    then has exactly two entries. Each of these entries is a dictionary of monovariate power
    sums as described for the parent class. These entries have respective keys 'upper' and
    'lower' for the first and second solver levels in the array received by self._update. See
    self._initialisePowerSums for details.

    * _powerSumsUpdater : callable
    Default: xmc.methodDefs_momentEstimator.powerSums.addPowerSumsAndCounter.
    Function which updates current power sums with additional power sums, and likewise for the
    sample count; called by self.update (see also its documentation). It must abide by the
    following specifications; see default function for an example.
    Input arguments:
    1. dictionary of power sums: DefaultDict[str, List[float]]
    2. current counter of samples: int.
    3. array of power sums: List[List[Union[List[float], int]]]
    4. keys of power sum (in argument 1) in the order of argument 3 : tuple (Default:
    ('1','2', ...))
    Output arguments:
    1. dictionary of power sums (updated): DefaultDict[str, List[float]]

    * _estimation: callable
    Default: xmc.methodDefs_momentEstimator.hStatistics.hStatisticsMonoPowerSums.
    See the description of this attribute of the parent class. The only difference is that the
    power sums are always monovariate, hence the different default value.

    """

    def __init__(self, **keywordArgs):
        # Parent constructor
        super().__init__(**keywordArgs)
        #
        # Methods
        powerSumsUpdaterName = keywordArgs.get(
            "powerSumsUpdater",
            "xmc.methodDefs_momentEstimator.powerSums.addPowerSumsAndCounter",
        )
        # Modules for methods to compute moments
        estimationName = keywordArgs.get(
            "estimationFunction",
            "xmc.methodDefs_momentEstimator.hStatistics.hStatisticsMonoPowerSums",
        )
        #
        # Non-empty initialisation
        # Enforce parallelisation
        if self._isUpdateParallel and powerSumsUpdaterName[-6:-1] != "_Task":
            powerSumsUpdaterName += "_Task"
        if self._isEstimationParallel and estimationName[-6:-1] != "_Task":
            estimationName += "_Task"
        self._powerSumsUpdater = dynamicImport(powerSumsUpdaterName)
        self._estimation = dynamicImport(estimationName)

    def update(self, samples):
        """
        Method updating the power sums and number of realizations, given new samples. The new
        data is passed through as a list containing power sums, up to a given order, and the
        number of samples from which it was computed. For example, let us discretise a scalar
        time signal u into M time steps to obtain a vector U of length M. A power sum of order
        a is computed as S_a = \sum_{i=1}^{M} U[i]^a.
        Samples will have the following shape:
        [[[[S1], [S2], M]]] ,
        where S1 and S2 are power sums of order one and two, respectively, and M is the total
        number of terms in the power sum (equivalent to number of time steps times realizations).

        Inputs:
        - self: an instance of the class
        - samples: nest list of dimension 4. Axes are: event, MC index, power, component.
        """
        # Synchronise if this estimator does support parallel computations
        if not self._isUpdateParallel:
            samples = get_value_from_remote(samples)
        # Deduce from self._indexSetDimension whether this is the initial update
        if self._indexSetDimension is None:
            # Guess missing information from sample set
            self._indexSetDimension = self._guessIndexSetDimension(samples[0])
            if self._variableDimension is None:
                self._variableDimension = self._guessVariableDimension(samples[0])
            if self.order is None:
                self.order = self._guessOrder(samples[0])
            # Initialise power sums and their updater if they are not set
            if self._powerSums is None or self._powerSumsUpdater is None:
                self._initialisePowerSums()
        #
        # Proceed to update
        self._powerSums, self._sampleCounter = self._powerSumsUpdater(
            self._powerSums, self._sampleCounter, samples
        )

    def _guessOrder(self, singleSample: List[List[float]]) -> int:
        """
        Guesses the statistical order (value of MultiCombinedMomentEstimator.order)
        from the power sums received.
        This is done under the assumption that the order of power sums must be twice
        the statistical order.
        """
        orderOfPowerSums = len(singleSample[0] - 1)
        return int(ceil(orderOfPowerSums / 2))

    def _initialisePowerSums(self):
        # Generate all integers between 0 and twice the maximal order of moments to compute.
        powerSumsKeys = tuple(str(power) for power in range(1, 2 * self.order + 1))
        # Generate list of dictionaries from these keys;
        # Initial value: NumPy array of zeros, in long format to avoid overflow.
        upperDict = {k: zeros(self._variableDimension, dtype=float64) for k in powerSumsKeys}
        if self._indexSetDimension == 0:
            self._powerSums = upperDict
        else:
            lowerDict = {
                k: zeros(self._variableDimension, dtype=float64) for k in powerSumsKeys
            }
            self._powerSums = {"upper": upperDict, "lower": lowerDict}
