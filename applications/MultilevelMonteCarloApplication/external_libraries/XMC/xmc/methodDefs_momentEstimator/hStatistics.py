from xmc.distributedEnvironmentFramework import ExaquteTask
from .types import PowerSumsDict, HStatistics, ListIndex
from numpy import ceil, log2


def hStatistics(
    component: ListIndex,
    powerSums: PowerSumsDict,
    numberOfSamples: int,
    order: int,
    isError: bool,
    isCentral: bool,
) -> HStatistics:
    """
    This is the generic function to compute h-statistics using the hStatistics module.

    It extracts the relevant power sums for the selected component(s), and calls the right
    function for the requested h-statistics. It abides by the specifications laid out in the
    documentation of xmc.momentEstimator.MultiMomentEstimator.
    """
    # Guess dimension d of the index set:
    # keys of dictionary of power sums have length 2^d
    indexSetDimension = int(ceil(log2(len(next(iter(powerSums))))))
    # Check that request can be fulfilled
    if indexSetDimension > 1 or order > 4:
        raise ValueError(
            f"Moments of order {order} are not supported "
            f"for a Monte Carlo index set of dimension {indexSetDimension}."
        )

    #
    # Prepare power sums
    if component != slice(None, None):
        # Extract power sums for chosen component
        powerSumsExcerpt = dict.fromkeys(powerSums.keys(), None)
        if hasattr(component, "__iter__"):
            # The variable component is an iterable of the indices of the components
            for key, val in powerSums.items():
                powerSumsExcerpt[key] = [val[i] for i in component]
        else:
            # The variable component is either an integer or a (strict) slice
            for key, val in powerSums.items():
                # If component is an integer, then powerSums[key] is a float
                # and the output of the function will be a signle float, and not a list.
                powerSumsExcerpt[key] = val[component]
        powerSums = powerSumsExcerpt
    # Assemble the name of function to compute the estimation
    functionName = (
        f"{'central' if isCentral else 'raw'}Moment"
        f"{'Error' if isError else ''}Order{order}Dimension{indexSetDimension}"
    )
    # Get this function
    statFun = globals().get(functionName)
    if not statFun:
        # The function does not exist in this module
        raise NotImplementedError(f"Function {functionName} not implemented")
    return statFun(powerSums, numberOfSamples).tolist()


@ExaquteTask()
def hStatistics_Task(
    component: ListIndex,
    powerSums: PowerSumsDict,
    numberOfSamples: int,
    order: int,
    isError: bool,
    isCentral: bool,
) -> HStatistics:
    """
    This is the parallelised version of hStatistic, from the same module. See its documentation.
    """
    return hStatistics(component, powerSums, numberOfSamples, order, isError, isCentral)


def hStatisticsMonoPowerSums(
    component: ListIndex,
    powerSums: PowerSumsDict,
    numberOfSamples: int,
    order: int,
    isError: bool,
    isCentral: bool,
) -> HStatistics:
    """
    Compute h-statistics from monovariate power sums. Obviously, this makes a difference only
    for differences of h-statistics (e.g. for MLMC). In this case, error estimation is not
    supported. This is useful for combined moment estimators, since they store only
    monovariate power sums.

    The input and output arguments are the same as for hStatistics, which it calls.
    """
    # Easy case first: single level
    if "lower" not in powerSums.keys():
        return hStatistics(component, powerSums, numberOfSamples, order, isError, isCentral)

    # The code below is for differences of h-statistics
    if isError:
        raise ValueError(
            "Error estimation is not implemented for differences of h-statistics "
            "computed from monovariate power sums."
        )

    hUpper = hStatistics(
        component, powerSums["upper"], numberOfSamples, order, isError, isCentral
    )
    hLower = hStatistics(
        component, powerSums["lower"], numberOfSamples, order, isError, isCentral
    )
    return hUpper - hLower


@ExaquteTask()
def hStatisticsMonoPowerSums_Task(
    component: ListIndex,
    powerSums: PowerSumsDict,
    numberOfSamples: int,
    order: int,
    isError: bool,
    isCentral: bool,
) -> HStatistics:
    """
    This is the parallelised version of hStatisticMonoPowerSums, from the same module.
    See its documentation.
    """
    return hStatisticsMonoPowerSums(
        component, powerSums, numberOfSamples, order, isError, isCentral
    )


# Organisation: Monte Carlo index dimension, moment order, {value, error}.

# Dimension 0

## Order 1

# TODO this is not a h-statistics. Is the module name misleading?
def rawMomentOrder1Dimension0(powerSums: PowerSumsDict, numberOfSamples: int) -> HStatistics:
    return powerSums["1"] / numberOfSamples


@ExaquteTask()
def rawMomentOrder1Dimension0_Task(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return rawMomentOrder1Dimension0(powerSums, numberOfSamples)


def rawMomentErrorOrder1Dimension0(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return centralMomentOrder2Dimension0(powerSums, numberOfSamples) / numberOfSamples


@ExaquteTask()
def rawMomentErrorOrder1Dimension0_Task(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return rawMomentErrorOrder1Dimension0(powerSums, numberOfSamples)


## Order 2


def centralMomentOrder2Dimension0(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return (numberOfSamples * powerSums["2"] - powerSums["1"] ** 2) / (
        (numberOfSamples - 1) * numberOfSamples
    )


@ExaquteTask()
def centralMomentOrder2Dimension0_Task(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return centralMomentOrder2Dimension0(powerSums, numberOfSamples)


def centralMomentErrorOrder2Dimension0(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    centralMoment2 = centralMomentOrder2Dimension0(powerSums, numberOfSamples)
    centralMoment4 = centralMomentOrder4Dimension0(powerSums, numberOfSamples)
    return centralMoment4 / numberOfSamples - (centralMoment2 ** 2) * (numberOfSamples - 3) / (
        numberOfSamples * (numberOfSamples - 1)
    )


@ExaquteTask()
def centralMomentErrorOrder2Dimension0_Task(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return centralMomentErrorOrder2Dimension0(powerSums, numberOfSamples)


## Order 3


def centralMomentOrder3Dimension0(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    # Source: MATHICSE technical report 23.2017, p. 5, 2nd unnumbered equation
    term1 = numberOfSamples ** 2 * powerSums["3"]
    term2 = 3.0 * numberOfSamples * powerSums["2"] * powerSums["1"]
    term3 = 2.0 * powerSums["1"] ** 3
    term4 = (numberOfSamples - 2.0) * (numberOfSamples - 1.0) * numberOfSamples
    return (term1 - term2 + term3) / term4


@ExaquteTask()
def centralMomentOrder3Dimension0_Task(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return centralMomentOrder3Dimension0(powerSums, numberOfSamples)


def centralMomentErrorOrder3Dimension0(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return _uglyD0O3E(powerSums, numberOfSamples)


@ExaquteTask()
def centralMomentErrorOrder3Dimension0_Task(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return centralMomentErrorOrder3Dimension0(powerSums, numberOfSamples)


## Order 4


def centralMomentOrder4Dimension0(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    # Source: MATHICSE technical report 23.2017, p. 5, third unnumbered equation
    term1 = (
        (-4.0 * numberOfSamples ** 2 + 8.0 * numberOfSamples - 12.0)
        * powerSums["3"]
        * powerSums["1"]
    )
    term2 = (
        numberOfSamples ** 3 - 2.0 * numberOfSamples ** 2 + 3.0 * numberOfSamples
    ) * powerSums["4"]
    term3 = 6.0 * numberOfSamples * powerSums["2"] * powerSums["1"] ** 2
    term4 = (9.0 - 6.0 * numberOfSamples) * powerSums["2"] ** 2
    term5 = 3.0 * powerSums["1"] ** 4
    term6 = (
        (numberOfSamples - 3.0)
        * (numberOfSamples - 2.0)
        * (numberOfSamples - 1.0)
        * numberOfSamples
    )
    return (term1 + term2 + term3 + term4 - term5) / term6


@ExaquteTask()
def centralMomentOrder4Dimension0_Task(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return centralMomentOrder4Dimension0(powerSums, numberOfSamples)


def centralMomentErrorOrder4Dimension0(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return _uglyD0O4E(powerSums, numberOfSamples)


@ExaquteTask()
def centralMomentErrorOrder4Dimension0_Task(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return centralMomentErrorOrder4Dimension0(powerSums, numberOfSamples)


# Dimension 1

# Order 1


def rawMomentOrder1Dimension1(powerSums: PowerSumsDict, numberOfSamples: int) -> HStatistics:
    # Select only power sums of differences
    return powerSums["01"] / numberOfSamples


@ExaquteTask()
def rawMomentOrder1Dimension1_Task(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return rawMomentOrder1Dimension1(powerSums, numberOfSamples)


def rawMomentErrorOrder1Dimension1(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    partialCentralMoment2 = (numberOfSamples * powerSums["02"] - powerSums["01"] ** 2) / (
        (numberOfSamples - 1) * numberOfSamples
    )
    # TODO should it not be powerSum01*powerSum10 instead of powerSum01**2?
    return partialCentralMoment2 / numberOfSamples


@ExaquteTask()
def rawMomentErrorOrder1Dimension1_Task(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return rawMomentErrorOrder1Dimension1(powerSums, numberOfSamples)


## Order 2


def centralMomentOrder2Dimension1(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    # Source: MATHICSE technical report 23.2017, p. 9, 2nd unnumbered equation above eq. (6a)
    return (numberOfSamples * powerSums["11"] - powerSums["10"] * powerSums["01"]) / (
        numberOfSamples * (numberOfSamples - 1)
    )


@ExaquteTask()
def centralMomentOrder2Dimension1_Task(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return centralMomentOrder2Dimension1(powerSums, numberOfSamples)


def centralMomentErrorOrder2Dimension1(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    # Source: MATHICSE technical report 23.2017, p. 12 eq. (7)
    term1 = numberOfSamples * (
        (-(numberOfSamples ** 2) + numberOfSamples + 2) * powerSums["11"] ** 2
        + (numberOfSamples * powerSums["22"] - 2 * powerSums["10"] * powerSums["12"])
        * (numberOfSamples - 1) ** 2
        + (numberOfSamples - 1) * powerSums["02"] * (powerSums["10"] ** 2 - powerSums["20"])
    )
    term2 = (powerSums["01"] ** 2) * (
        (6 - 4 * numberOfSamples) * powerSums["10"] ** 2
        + (numberOfSamples - 1) * numberOfSamples * powerSums["20"]
    )
    term3 = (
        2
        * numberOfSamples
        * powerSums["01"]
        * (
            powerSums["21"] * (numberOfSamples - 1) ** 2
            + (5 - 3 * numberOfSamples) * powerSums["10"] * powerSums["11"]
        )
    )
    term4 = (
        (numberOfSamples - 3)
        * (numberOfSamples - 2)
        * (numberOfSamples - 1) ** 2
        * numberOfSamples ** 2
    )
    return (term1 + term2 - term3) / term4


@ExaquteTask()
def centralMomentErrorOrder2Dimension1_Task(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return centralMomentErrorOrder2Dimension1(powerSums, numberOfSamples)


## Order 3


def centralMomentOrder3Dimension1(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    # Source: MATHICSE technical report 23.2017, p. 9, unnumbered equation above (6a)
    term1 = numberOfSamples ** 2 * (powerSums["03"] + 3 * powerSums["21"])
    term2 = 3.0 * numberOfSamples * powerSums["01"] * (powerSums["02"] + powerSums["20"])
    term3 = 6.0 * numberOfSamples * powerSums["10"] * powerSums["11"]
    term4 = 2.0 * powerSums["01"] ** 3
    term5 = 6.0 * numberOfSamples * powerSums["01"] * powerSums["10"] ** 2
    term6 = 4.0 * (numberOfSamples - 2.0) * (numberOfSamples - 1.0) * numberOfSamples
    return (-term1 + term2 + term3 - term4 - term5) / term6


@ExaquteTask()
def centralMomentOrder3Dimension1_Task(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return centralMomentOrder3Dimension1(powerSums, numberOfSamples)


def centralMomentErrorOrder3Dimension1(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return _uglyD1O3E(powerSums, numberOfSamples)


@ExaquteTask()
def centralMomentErrorOrder3Dimension1_Task(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return centralMomentErrorOrder3Dimension1(powerSums, numberOfSamples)


## Order 4


def centralMomentOrder4Dimension1(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return _uglyD1O4V(powerSums, numberOfSamples)


@ExaquteTask()
def centralMomentOrder4Dimension1_Task(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return centralMomentOrder4Dimension1(powerSums, numberOfSamples)


def centralMomentErrorOrder4Dimension1(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return _uglyD1O4E(powerSums, numberOfSamples)


@ExaquteTask()
def centralMomentErrorOrder4Dimension1_Task(
    powerSums: PowerSumsDict, numberOfSamples: int
) -> HStatistics:
    return centralMomentErrorOrder4Dimension1(powerSums, numberOfSamples)


###############################################################################################
#                                  Machine-generated code                                     #
###############################################################################################

# Dimension 0

## Order 3


def _uglyD0O3E(powerSums: PowerSumsDict, numberOfSamples: int) -> HStatistics:
    # Dimension 0, order 3, error
    # Source: MATHICSE technical report 23.2017, table 1, p. 7
    t1 = numberOfSamples ** 2
    t2 = t1 * numberOfSamples
    t3 = t1 ** 2
    t6 = t3 * t1
    t13 = powerSums["3"] ** 2
    t15 = t3 * numberOfSamples
    t16 = powerSums["1"] ** 2
    t20 = t15 * powerSums["1"]
    t21 = powerSums["2"] * powerSums["3"]
    t24 = powerSums["2"] ** 2
    t25 = t24 * powerSums["2"]
    t28 = t16 * powerSums["1"]
    t32 = t3 * t16
    t35 = t16 ** 2
    t39 = t35 * t16
    t44 = (
        21 * t15 * t16 * powerSums["4"]
        + t3 * t2 * powerSums["6"]
        + 108 * t2 * t35 * powerSums["2"]
        - 48 * t3 * t28 * powerSums["3"]
        - 6 * t6 * powerSums["1"] * powerSums["5"]
        - 6 * t6 * powerSums["2"] * powerSums["4"]
        - 36 * t1 * t39
        - t6 * t13
        + 9 * t15 * t25
        + 30 * t20 * t21
        - 72 * t32 * t24
        - 5 * t6 * powerSums["6"]
    )
    t54 = t3 * powerSums["1"]
    t62 = t2 * t16
    t77 = (
        -540 * t1 * t35 * powerSums["2"]
        + 33 * t15 * powerSums["2"] * powerSums["4"]
        + 216 * t2 * t28 * powerSums["3"]
        - 42 * t3 * powerSums["2"] * powerSums["4"]
        + 180 * numberOfSamples * t39
        - 4 * t15 * t13
        + 13 * t15 * powerSums["6"]
        + 30 * t20 * powerSums["5"]
        - 108 * t54 * t21
        + 378 * t62 * t24
        - 72 * t3 * t25
        - 108 * t32 * powerSums["4"]
        - 78 * t54 * powerSums["5"]
    )
    t83 = t2 * powerSums["1"]
    t91 = t1 * t16
    t109 = (
        720 * numberOfSamples * t35 * powerSums["2"]
        - 264 * t1 * t28 * powerSums["3"]
        - 75 * t2 * powerSums["2"] * powerSums["4"]
        - 40 * t2 * t13
        + 41 * t3 * t13
        + 213 * t2 * t25
        - 78 * t83 * t21
        - 522 * t91 * t24
        - 23 * t3 * powerSums["6"]
        + 237 * t62 * powerSums["4"]
        + 138 * t83 * powerSums["5"]
        - 270 * t91 * powerSums["4"]
        - 240 * t39
    )
    t110 = t1 * powerSums["1"]
    t127 = numberOfSamples * powerSums["1"]
    t141 = (
        120 * numberOfSamples * t16 * powerSums["4"]
        - 120 * numberOfSamples * powerSums["2"] * powerSums["4"]
        + 210 * t1 * powerSums["2"] * powerSums["4"]
        + 80 * numberOfSamples * t13
        + 120 * numberOfSamples * t25
        - 100 * t1 * t13
        - 270 * t1 * t25
        - 8 * t1 * powerSums["6"]
        + 540 * t110 * t21
        - 132 * t110 * powerSums["5"]
        - 240 * t127 * t21
        + 48 * t127 * powerSums["5"]
        + 22 * t2 * powerSums["6"]
    )
    t155 = (-2 + numberOfSamples) ** 2
    t158 = (-1 + numberOfSamples) ** 2
    uglyD0O3E = (
        (t44 + t77 + t109 + t141)
        / numberOfSamples
        / (-5 + numberOfSamples)
        / (-4 + numberOfSamples)
        / (-3 + numberOfSamples)
        / t155
        / t158
    )
    return uglyD0O3E


## Order 4


def _uglyD0O4E(powerSums: PowerSumsDict, numberOfSamples: int) -> HStatistics:
    # Dimension 0, order 4, error
    # Source: MATHICSE technical report 23.2017, appendix A, p. 27
    t1 = numberOfSamples ** 2
    t2 = t1 ** 2
    t3 = t2 ** 2
    t6 = t3 * numberOfSamples
    t13 = powerSums["4"] ** 2
    t15 = powerSums["1"] ** 2
    t19 = t3 * powerSums["1"]
    t20 = powerSums["2"] * powerSums["5"]
    t23 = powerSums["3"] * powerSums["4"]
    t26 = t3 * powerSums["2"]
    t27 = powerSums["3"] ** 2
    t30 = t1 * numberOfSamples
    t31 = t2 * t30
    t32 = t15 * powerSums["1"]
    t36 = t31 * t15
    t37 = powerSums["2"] * powerSums["4"]
    t42 = t31 * powerSums["1"]
    t43 = powerSums["2"] ** 2
    t44 = t43 * powerSums["3"]
    t47 = t2 * t1
    t48 = t15 ** 2
    t52 = t47 * t32
    t53 = powerSums["2"] * powerSums["3"]
    t56 = t47 * t15
    t57 = t43 * powerSums["2"]
    t60 = t2 * numberOfSamples
    t61 = t48 * powerSums["1"]
    t65 = t60 * t48
    t68 = t48 * t15
    t72 = t48 ** 2
    t75 = (
        t3 * t1 * powerSums["8"]
        - 8 * t6 * powerSums["1"] * powerSums["7"]
        - 8 * t6 * powerSums["3"] * powerSums["5"]
        - t6 * t13
        + 28 * t3 * t15 * powerSums["6"]
        + 24 * t19 * t20
        + 48 * t19 * t23
        + 16 * t26 * t27
        - 72 * t31 * t32 * powerSums["5"]
        - 132 * t36 * t37
        - 112 * t36 * t27
        - 96 * t42 * t44
        + 156 * t47 * t48 * powerSums["4"]
        + 528 * t52 * t53
        + 144 * t56 * t57
        - 336 * t60 * t61 * powerSums["3"]
        - 612 * t65 * t43
        + 576 * t2 * t68 * powerSums["2"]
        - 144 * t30 * t72
    )
    t95 = t31 * powerSums["2"]
    t104 = t47 * powerSums["1"]
    t107 = t43 ** 2
    t112 = t60 * t32
    t115 = t60 * t15
    t121 = t2 * t48
    t124 = (
        -17 * t6 * powerSums["8"]
        + 136 * t19 * powerSums["7"]
        - 52 * t26 * powerSums["6"]
        + 160 * t3 * powerSums["3"] * powerSums["5"]
        + t3 * t13
        - 424 * t36 * powerSums["6"]
        - 168 * t42 * t20
        - 808 * t42 * t23
        + 132 * t31 * t43 * powerSums["4"]
        - 208 * t95 * t27
        + 960 * t52 * powerSums["5"]
        + 1368 * t56 * t37
        + 1824 * t56 * t27
        + 720 * t104 * t44
        - 72 * t47 * t107
        - 1884 * t65 * powerSums["4"]
        - 6432 * t112 * t53
        - 792 * t115 * t57
        + 4080 * t2 * t61 * powerSums["3"]
        + 6012 * t121 * t43
    )
    t151 = t47 * powerSums["2"]
    t160 = t60 * powerSums["1"]
    t167 = t2 * t32
    t170 = t2 * t15
    t173 = (
        -6048 * t30 * t68 * powerSums["2"]
        + 1512 * t1 * t72
        + 140 * t3 * powerSums["8"]
        - 1120 * t42 * powerSums["7"]
        + 808 * t95 * powerSums["6"]
        - 1504 * t31 * powerSums["3"] * powerSums["5"]
        + 148 * t31 * t13
        + 3112 * t56 * powerSums["6"]
        - 336 * t104 * t20
        + 6336 * t104 * t23
        - 1620 * t47 * t43 * powerSums["4"]
        + 880 * t151 * t27
        - 6000 * t112 * powerSums["5"]
        - 5424 * t115 * t37
        - 13552 * t115 * t27
        + 1200 * t160 * t44
        + 612 * t60 * t107
        + 10212 * t121 * powerSums["4"]
        + 32736 * t167 * t53
        - 4248 * t170 * t57
    )
    t177 = t30 * t48
    t205 = t60 * powerSums["2"]
    t214 = t2 * powerSums["1"]
    t221 = (
        -21264 * t30 * t61 * powerSums["3"]
        - 18180 * t177 * t43
        + 22752 * t1 * t68 * powerSums["2"]
        - 5688 * numberOfSamples * t72
        - 650 * t31 * powerSums["8"]
        + 5200 * t104 * powerSums["7"]
        - 5656 * t151 * powerSums["6"]
        + 8512 * t47 * powerSums["3"] * powerSums["5"]
        - 1694 * t47 * t13
        - 12544 * t115 * powerSums["6"]
        + 8400 * t160 * t20
        - 29008 * t160 * t23
        + 6972 * t60 * t43 * powerSums["4"]
        + 368 * t205 * t27
        + 19488 * t167 * powerSums["5"]
        + 8568 * t170 * t37
        + 57648 * t170 * t27
        - 30096 * t214 * t44
        + 828 * t2 * t107
        - 28644 * t177 * powerSums["4"]
    )
    t224 = t30 * t32
    t227 = t30 * t15
    t233 = t1 * t48
    t260 = t2 * powerSums["2"]
    t269 = (
        -30240 * numberOfSamples * t68 * powerSums["2"]
        + 57552 * t1 * t61 * powerSums["3"]
        - 10908 * t2 * t43 * powerSums["4"]
        - 32296 * t60 * powerSums["3"] * powerSums["5"]
        + 9715 * t60 * t13
        - 15016 * t160 * powerSums["7"]
        + 29548 * t170 * powerSums["6"]
        - 41160 * t214 * t20
        + 23008 * t205 * powerSums["6"]
        + 83760 * t214 * t23
        - 86592 * t224 * t53
        - 31656 * t224 * powerSums["5"]
        - 152512 * t227 * t27
        - 924 * t227 * t37
        + 41832 * t227 * t57
        + 2196 * t233 * t43
        - 15008 * t260 * t27
        + 1877 * t47 * powerSums["8"]
        + 7560 * t72
    )
    t270 = t30 * powerSums["1"]
    t277 = t1 * t32
    t280 = t1 * t15
    t286 = numberOfSamples * t48
    t309 = t30 * powerSums["2"]
    t318 = (
        -83232 * numberOfSamples * t61 * powerSums["3"]
        + 84640 * t2 * powerSums["3"] * powerSums["5"]
        + 240 * t30 * t43 * powerSums["4"]
        - 23724 * t30 * t107
        - 33983 * t2 * t13
        + 94584 * t270 * t20
        + 28936 * t214 * powerSums["7"]
        - 43192 * t227 * powerSums["6"]
        - 151336 * t270 * t23
        + 40032 * t233 * powerSums["4"]
        - 58084 * t260 * powerSums["6"]
        + 253872 * t280 * t27
        + 48800 * t309 * t27
        + 133680 * t270 * t44
        + 128016 * t277 * t53
        + 23328 * t277 * powerSums["5"]
        - 9936 * t280 * t37
        - 105624 * t280 * t57
        + 62424 * t286 * t43
        - 3617 * t60 * powerSums["8"]
    )
    t320 = t1 * powerSums["1"]
    t327 = numberOfSamples * t32
    t330 = numberOfSamples * t15
    t357 = t1 * powerSums["2"]
    t366 = (
        6480 * t1 * t43 * powerSums["4"]
        - 144176 * t30 * powerSums["3"] * powerSums["5"]
        + 89748 * t1 * t107
        + 70850 * t30 * t13
        + 4534 * t2 * powerSums["8"]
        - 103680 * t320 * t20
        + 154080 * t320 * t23
        - 245088 * t330 * t27
        - 63072 * t357 * t27
        - 36272 * t270 * powerSums["7"]
        + 37584 * t280 * powerSums["6"]
        - 24192 * t286 * powerSums["4"]
        + 89368 * t309 * powerSums["6"]
        - 293760 * t320 * t44
        - 102816 * t327 * t53
        - 6048 * t327 * powerSums["5"]
        + 15120 * t330 * t37
        + 81648 * t330 * t57
        - 45360 * t48 * t43
        + 60480 * t61 * powerSums["3"]
    )
    t367 = numberOfSamples * powerSums["1"]
    t392 = numberOfSamples * powerSums["2"]
    t412 = (
        3024 * numberOfSamples * t43 * powerSums["4"]
        - 52416 * numberOfSamples * powerSums["3"] * powerSums["5"]
        + 137088 * t1 * powerSums["3"] * powerSums["5"]
        - 181440 * powerSums["1"] * t43 * powerSums["3"]
        - 132192 * numberOfSamples * t107
        + 32760 * numberOfSamples * t13
        - 76356 * t1 * t13
        + 936 * t1 * powerSums["8"]
        + 120960 * t15 * t27
        + 42336 * t367 * t20
        - 74592 * t367 * t23
        + 28224 * t392 * t27
        - 3204 * t30 * powerSums["8"]
        + 25632 * t320 * powerSums["7"]
        - 14112 * t330 * powerSums["6"]
        - 75600 * t357 * powerSums["6"]
        + 352512 * t367 * t44
        - 7488 * t367 * powerSums["7"]
        + 26208 * t392 * powerSums["6"]
        + 68040 * t107
    )
    t430 = (-3 + numberOfSamples) ** 2
    t433 = (-2 + numberOfSamples) ** 2
    t437 = (-1 + numberOfSamples) ** 2
    uglyD0O4E = (
        (t75 + t124 + t173 + t221 + t269 + t318 + t366 + t412)
        / numberOfSamples
        / (-7 + numberOfSamples)
        / (-6 + numberOfSamples)
        / (-5 + numberOfSamples)
        / (-4 + numberOfSamples)
        / t430
        / t433
        / t437
    )
    return uglyD0O4E


# Dimension 1

## Order 3


def _uglyD1O3E(powerSums: PowerSumsDict, numberOfSamples: int) -> HStatistics:
    # Dimension 1, order 3, error
    # Source: MATHICSE technical report 23.2017, appendix B, p. 28
    t1 = numberOfSamples ** 2
    t2 = 3.0 * t1
    t4 = t2 - 15 * numberOfSamples + 20.0
    t5 = powerSums["01"] ** 2
    t6 = t5 ** 2
    t12 = powerSums["10"] ** 2
    t15 = -3.0 + numberOfSamples
    t16 = t15 ** 2
    t27 = 2.0 * t1 - 9.0 * numberOfSamples + 11.0
    t33 = numberOfSamples * t1
    t34 = 5.0 * t1
    t35 = 6.0 * numberOfSamples
    t36 = t33 - t34 + t35 + 2.0
    t46 = powerSums["02"] ** 2
    t49 = t1 ** 2
    t53 = 90.0 * numberOfSamples
    t58 = 8.0 * t1
    t62 = 3.0 * t33
    t64 = 30.0 * numberOfSamples
    t70 = t12 ** 2
    t79 = -2.0 + numberOfSamples
    t80 = -1.0 + numberOfSamples
    t81 = t80 ** 2
    t82 = t79 * t81
    t83 = numberOfSamples * powerSums["10"]
    t89 = powerSums["11"] ** 2
    t92 = 2.0 * t33
    t94 = 7.0 * numberOfSamples
    t96 = powerSums["20"] ** 2
    t100 = t1 - 3.0 * numberOfSamples + 2.0
    t115 = numberOfSamples * t49
    t122 = powerSums["10"] * t12
    t144 = (
        88.0 * numberOfSamples * powerSums["23"]
        + 66.0 * numberOfSamples * powerSums["41"]
        - 92.0 * t1 * powerSums["23"]
        + 3.0 * t115 * powerSums["41"]
        + 24.0 * powerSums["03"] * t12
        + 48.0 * t12 * powerSums["21"]
        - 696.0 * t122 * powerSums["11"]
        + 52.0 * t33 * powerSums["23"]
        + 39.0 * t33 * powerSums["41"]
        + 16.0 * powerSums["03"] * powerSums["20"]
        - 64.0 * powerSums["10"] * powerSums["13"]
        - 96.0 * powerSums["10"] * powerSums["31"]
        + 96.0 * powerSums["11"] * powerSums["12"]
        - 32.0 * powerSums["23"]
    )
    t165 = 3.0 * t49
    t174 = t49 * powerSums["10"]
    t177 = t33 * powerSums["10"]
    t180 = t1 * powerSums["10"]
    t185 = t49 * powerSums["11"]
    t188 = t1 * powerSums["11"]
    t191 = t33 * powerSums["11"]
    t194 = (
        -69.0 * t1 * powerSums["41"]
        + 72.0 * powerSums["20"] * powerSums["21"]
        - 20.0 * t49 * powerSums["23"]
        + 4.0 * t115 * powerSums["23"]
        + 48.0 * powerSums["11"] * powerSums["30"]
        - 15.0 * t49 * powerSums["41"]
        + powerSums["02"]
        * (
            (-5.0 * t49 + 18.0 * t33 + 13.0 * t1 - t53 + 40.0) * powerSums["03"]
            + 36.0 * t36 * powerSums["10"] * powerSums["11"]
            - 3.0 * (t165 - 14.0 * t33 + t1 + 50.0 * numberOfSamples - 16.0) * powerSums["21"]
        )
        - 24.0 * powerSums["41"]
        - 18.0 * t174 * powerSums["31"]
        + 96.0 * t177 * powerSums["31"]
        - 210.0 * t180 * powerSums["31"]
        + 228.0 * t83 * powerSums["31"]
        - 6.0 * t185 * powerSums["30"]
        - 6.0 * t188 * powerSums["30"]
        + 24.0 * t191 * powerSums["30"]
    )
    t196 = numberOfSamples * powerSums["11"]
    t201 = t1 * t122
    t207 = t33 * powerSums["03"]
    t210 = numberOfSamples * powerSums["03"]
    t217 = numberOfSamples * t122
    t220 = t1 * powerSums["03"]
    t231 = (
        -60.0 * t196 * powerSums["30"]
        + 90.0 * t188 * powerSums["12"]
        - 96.0 * t201 * powerSums["11"]
        + t81 * (t33 - t2 + t35 - 8.0) * powerSums["05"]
        + 12.0 * t207 * t12
        + 72.0 * t210 * t12
        - 180.0 * t196 * powerSums["12"]
        - 6.0 * t185 * powerSums["12"]
        + 504.0 * t217 * powerSums["11"]
        - 60.0 * t220 * t12
        + 132 * t83 * powerSums["13"]
        - 106.0 * t180 * powerSums["13"]
        - 10.0 * t174 * powerSums["13"]
        + 48.0 * t177 * powerSums["13"]
    )
    t235 = t49 * powerSums["20"]
    t238 = t33 * t12
    t241 = powerSums["11"] * powerSums["20"]
    t246 = t33 * powerSums["20"]
    t249 = t1 * t12
    t256 = t1 * powerSums["20"]
    t259 = numberOfSamples * t12
    t266 = numberOfSamples * powerSums["20"]
    t269 = powerSums["10"] * powerSums["11"]
    t272 = (
        -5.0 * t49 * powerSums["03"] * powerSums["20"]
        + 60.0 * t177 * t241
        - 324.0 * t180 * t241
        + 30.0 * t207 * powerSums["20"]
        - 30.0 * t210 * powerSums["20"]
        - 35.0 * t220 * powerSums["20"]
        - 9.0 * t235 * powerSums["21"]
        + 48.0 * t238 * powerSums["21"]
        + 480.0 * t83 * t241
        + 30.0 * t246 * powerSums["21"]
        - 228.0 * t249 * powerSums["21"]
        + 45.0 * t256 * powerSums["21"]
        + 276.0 * t259 * powerSums["21"]
        - 210.0 * t266 * powerSums["21"]
        - 72.0 * t269 * powerSums["20"]
    )
    t277 = powerSums["21"] ** 2
    t279 = powerSums["12"] ** 2
    t286 = powerSums["03"] ** 2
    t323 = (
        432.0 * t277
        + 288.0 * t279
        - (t115 + 4.0 * t49 - 41.0 * t33 + 40.0 * t1 + 100.0 * numberOfSamples - 80.0) * t286
        - 8.0 * numberOfSamples * powerSums["06"]
        + 22.0 * t1 * powerSums["06"]
        - 23.0 * t33 * powerSums["06"]
        + 864.0 * t89 * powerSums["20"]
        - 138.0 * t33 * powerSums["24"]
        - 432.0 * powerSums["20"] * powerSums["22"]
        + 6.0
        * powerSums["03"]
        * (
            4.0 * (t49 - t62 - t58 + t64 - 8.0) * powerSums["10"] * powerSums["11"]
            - (t115 - 2.0 * t49 - 17.0 * t33 + 34.0 * t1 + 40.0 * numberOfSamples - 32.0)
            * powerSums["21"]
        )
        + 3.0 * (t165 - 24.0 * t33 + 71.0 * t1 - t53 + 40.0) * powerSums["02"] * t46
        - 540.0 * numberOfSamples * t277
        - 48.0 * powerSums["04"] * powerSums["20"]
    )
    t324 = t1 * t49
    t351 = (
        -72.0 * numberOfSamples * powerSums["42"]
        + 132.0 * t1 * powerSums["24"]
        - 30.0 * t115 * powerSums["24"]
        - 45.0 * t115 * powerSums["42"]
        - 144.0 * t12 * t89
        + 432.0 * t12 * powerSums["22"]
        + 144.0 * t122 * powerSums["12"]
        + 6.0 * t324 * powerSums["24"]
        + 9.0 * t324 * powerSums["42"]
        + 78.0 * t49 * powerSums["24"]
        + 96.0 * powerSums["10"] * powerSums["14"]
        - 384.0 * powerSums["11"] * powerSums["13"]
        + 288.0 * powerSums["12"] * powerSums["30"]
    )
    t356 = powerSums["11"] * powerSums["21"]
    t359 = powerSums["12"] * powerSums["20"]
    t381 = (
        72.0 * t49 * t89 * powerSums["20"]
        - 9.0 * t115 * t277
        - 5.0 * t115 * powerSums["06"]
        + 144.0 * t174 * t356
        + 36.0 * t174 * t359
        - 648.0 * t177 * t356
        - 72.0 * t177 * t359
        - 12.0 * t191 * powerSums["13"]
        + 144.0 * t33 * t279
        - 36.0 * t49 * t279
        + t324 * powerSums["06"]
        + 13.0 * t49 * powerSums["06"]
        + 117.0 * t49 * powerSums["42"]
    )
    t410 = (
        -360.0 * numberOfSamples * t279
        - 48.0 * numberOfSamples * powerSums["24"]
        - 324.0 * t1 * t277
        - 36.0 * t1 * t279
        + 198.0 * t1 * powerSums["42"]
        + 48.0 * powerSums["04"] * t12
        + 72.0 * t180 * t356
        - 252.0 * t180 * t359
        + 225.0 * t33 * t277
        - 207.0 * t33 * powerSums["42"]
        + 2160.0 * t83 * t356
        + 720.0 * t83 * t359
        + 288.0 * powerSums["10"] * powerSums["32"]
        - 576.0 * powerSums["11"] * powerSums["31"]
    )
    t442 = (
        -360.0 * numberOfSamples * powerSums["12"] * powerSums["30"]
        - 36.0 * t49 * powerSums["12"] * powerSums["30"]
        - 432.0 * powerSums["10"] * powerSums["12"] * powerSums["20"]
        + 48.0 * t185 * powerSums["13"]
        + 216.0 * t185 * powerSums["31"]
        - 312.0 * t188 * powerSums["13"]
        + 672.0 * t196 * powerSums["13"]
        + 288.0 * t201 * powerSums["12"]
        - 216.0 * t238 * t89
        - 360.0 * t238 * powerSums["22"]
        - 90.0 * t246 * powerSums["22"]
        - 306.0 * t256 * powerSums["22"]
        - 1728.0 * t259 * t89
    )
    t458 = 5.0 * numberOfSamples
    t500 = (
        (t92 - t34 - t458 + 20.0) * powerSums["04"]
        - 12.0 * t70
        - 12.0 * (t1 - t458 + 8.0) * t89
        + 12.0 * (1.0 + numberOfSamples) * t12 * powerSums["20"]
        - 36.0 * t96
        + 15.0 * numberOfSamples * t96
        - 3.0 * t1 * t96
        + 48.0 * powerSums["22"]
        - 12.0 * numberOfSamples * powerSums["22"]
        - 18.0 * t1 * powerSums["22"]
        + 6.0 * t33 * powerSums["22"]
        + 3.0
        * powerSums["10"]
        * (-4.0 * (t1 - numberOfSamples - 4.0) * powerSums["12"] - 8.0 * t80 * powerSums["30"])
        + 12.0 * powerSums["40"]
        - 3.0 * numberOfSamples * powerSums["40"]
        + 3.0 * t1 * powerSums["40"]
    )
    t507 = t1 * powerSums["04"]
    t512 = t33 * powerSums["04"]
    t515 = t49 * powerSums["04"]
    t518 = t115 * powerSums["11"]
    t521 = (
        756.0 * t266 * powerSums["22"]
        + 144.0 * t33 * powerSums["12"] * powerSums["30"]
        - 36.0 * t1 * powerSums["12"] * powerSums["30"]
        + 90.0 * t235 * powerSums["22"]
        - 324.0 * t191 * powerSums["31"]
        + 18.0 * t100 * t46 * (-2.0 * t15 * t12 + (t1 - t458 + 4.0) * powerSums["20"])
        - 3.0 * t100 * powerSums["02"] * t500
        + 276.0 * t180 * powerSums["14"]
        - 360.0 * t217 * powerSums["12"]
        + 156.0 * t507 * t12
        + 1224.0 * t249 * t89
        - 72.0 * t512 * t12
        + 12.0 * t515 * t12
        - 36.0 * t518 * powerSums["31"]
    )
    t526 = t115 * powerSums["10"]
    t555 = (
        -6.0 * t115 * powerSums["04"] * powerSums["20"]
        - 18.0 * t115 * powerSums["20"] * powerSums["22"]
        + 72.0 * t49 * t12 * powerSums["22"]
        - 72.0 * t33 * t122 * powerSums["12"]
        - 576.0 * t33 * t89 * powerSums["20"]
        + 60.0 * t174 * powerSums["14"]
        + 180.0 * t174 * powerSums["32"]
        - 156.0 * t177 * powerSums["14"]
        - 78.0 * t512 * powerSums["20"]
        + 42.0 * t515 * powerSums["20"]
        - 12.0 * t518 * powerSums["13"]
        - 12.0 * t526 * powerSums["14"]
        - 36.0 * t526 * powerSums["32"]
    )
    t569 = numberOfSamples * powerSums["04"]
    t587 = (
        -2016.0 * numberOfSamples * t89 * powerSums["20"]
        + 1656.0 * t1 * t89 * powerSums["20"]
        - 144.0 * t569 * t12
        - 468.0 * t177 * powerSums["32"]
        + 828.0 * t180 * powerSums["32"]
        - 288.0 * t188 * powerSums["31"]
        + 1008.0 * t196 * powerSums["31"]
        + 792.0 * t249 * powerSums["22"]
        - 936.0 * t259 * powerSums["22"]
        - 864.0 * t269 * powerSums["21"]
        + 6.0 * t507 * powerSums["20"]
        + 84.0 * t569 * powerSums["20"]
        - 264.0 * t83 * powerSums["14"]
        - 792.0 * t83 * powerSums["32"]
    )
    t600 = t79 ** 2
    uglyD1O3E = (
        (
            -12.0 * t4 * t5 * t6
            + 36.0
            * t6
            * (
                2.0 * t16 * numberOfSamples * powerSums["20"]
                + numberOfSamples * t4 * powerSums["02"]
                - 2.0 * t4 * t12
            )
            - 24.0
            * numberOfSamples
            * powerSums["01"]
            * t5
            * (
                numberOfSamples * t27 * powerSums["03"]
                - 6.0 * t27 * powerSums["10"] * powerSums["11"]
                + 3.0 * t36 * powerSums["21"]
            )
            + 3.0
            * t5
            * (
                -6.0 * t1 * (4.0 * t1 - 21.0 * numberOfSamples + 29.0) * t46
                + numberOfSamples
                * (7.0 * t49 - 36.0 * t33 + 79.0 * t1 - t53 + 40.0)
                * powerSums["04"]
                - 12.0
                * numberOfSamples
                * powerSums["02"]
                * (
                    (-t58 + 42.0 * numberOfSamples - 58.0) * t12
                    + (t62 - 19.0 * t1 + t64 - 2.0) * powerSums["20"]
                )
                - 36.0 * t4 * t70
                + 24.0
                * numberOfSamples
                * (t34 - 24.0 * numberOfSamples + 31.0)
                * t12
                * powerSums["20"]
                - 24.0 * t82 * t83 * (2.0 * powerSums["12"] + powerSums["30"])
                + 3.0
                * numberOfSamples
                * (
                    -8.0 * t82 * t89
                    - 2.0 * (t92 - 9.0 * t1 + t94 + 12.0) * t96
                    + t100
                    * (
                        2.0 * (t2 - t94 + 8.0) * powerSums["22"]
                        + (t1 - numberOfSamples + 4.0) * powerSums["40"]
                    )
                )
            )
            - 6.0 * numberOfSamples * powerSums["01"] * (t144 + t194 + t231 + t272)
            + numberOfSamples * (t323 + t351 + t381 + t410 + t442 + t521 + t555 + t587)
        )
        / (-5.0 + numberOfSamples)
        / (-4.0 + numberOfSamples)
        / t15
        / t600
        / t81
        / numberOfSamples
        / 16.0
    )
    return uglyD1O3E


## Order 4


def _uglyD1O4V(powerSums: PowerSumsDict, numberOfSamples: int) -> HStatistics:
    # Dimension 1, order 4, value
    # Source: trust me
    t1 = (
        2
        * numberOfSamples
        * (numberOfSamples - 1)
        * (numberOfSamples - 2)
        * (numberOfSamples - 3)
    )
    t2 = numberOfSamples ** 2 - 2 * numberOfSamples + 3
    t3 = 9 - 6 * numberOfSamples
    t4 = -3 * t2 * powerSums["01"] * powerSums["12"]
    t5 = -3 * t2 * powerSums["10"] * powerSums["21"]
    t6 = numberOfSamples * t2 * powerSums["13"]
    t7 = numberOfSamples * t2 * powerSums["31"]
    t8 = t3 * powerSums["11"] * powerSums["20"]
    t9 = t3 * powerSums["11"] * powerSums["02"]
    t10 = (
        3
        * numberOfSamples
        * powerSums["10"]
        * powerSums["01"]
        * (powerSums["20"] + powerSums["02"])
    )
    t11 = -t2 * powerSums["01"] * powerSums["30"]
    t12 = -t2 * powerSums["10"] * powerSums["03"]
    t13 = (
        3
        * (powerSums["01"] ** 2 + powerSums["10"] ** 2)
        * (numberOfSamples * powerSums["11"] - powerSums["10"] * powerSums["01"])
    )
    uglyD1O4V = (t4 + t5 + t6 + t7 + t8 + t9 + t10 + t11 + t12 + t13) / t1
    return uglyD1O4V


def _uglyD1O4E(powerSums: PowerSumsDict, numberOfSamples: int) -> HStatistics:
    # Dimension 1, order 4, error
    # Source: cmlmc-py library (unreleased)
    t1 = numberOfSamples ** 2
    t2 = t1 ** 2
    t3 = t2 ** 2
    t4 = t3 * powerSums["10"]
    t5 = powerSums["14"] * powerSums["20"]
    t8 = t1 * numberOfSamples
    t9 = t2 * t8
    t10 = powerSums["10"] ** 2
    t11 = t9 * t10
    t12 = powerSums["20"] * powerSums["22"]
    t15 = t2 * t1
    t16 = t10 * powerSums["10"]
    t17 = t15 * t16
    t18 = powerSums["11"] * powerSums["21"]
    t21 = powerSums["12"] * powerSums["20"]
    t24 = t9 * powerSums["10"]
    t25 = powerSums["12"] * powerSums["40"]
    t30 = t15 * t10
    t33 = t15 * powerSums["10"]
    t34 = powerSums["20"] ** 2
    t35 = powerSums["12"] * t34
    t38 = t2 * numberOfSamples
    t39 = t38 * t16
    t49 = t38 * t10
    t52 = t38 * powerSums["10"]
    t55 = t2 * t16
    t64 = t2 * t10
    t67 = t2 * powerSums["10"]
    t70 = t8 * t16
    t73 = (
        -198 * t49 * t12
        - 666 * t64 * t12
        + 12258 * t55 * t18
        - 34074 * t70 * t18
        - 774 * t55 * t21
        - 72 * t33 * t25
        + 12 * t52 * t25
        + 312 * t33 * t5
        + 1008 * t52 * t35
        - 4284 * t67 * t35
        - 252 * t52 * t5
    )
    t81 = t8 * t10
    t84 = t8 * powerSums["10"]
    t87 = t1 * t16
    t96 = t1 * t10
    t100 = t1 * powerSums["10"]
    t103 = numberOfSamples * t16
    t106 = numberOfSamples * t10
    t107 = powerSums["11"] ** 2
    t108 = t107 * powerSums["20"]
    t115 = powerSums["20"] * powerSums["32"]
    t120 = numberOfSamples * powerSums["10"]
    t129 = (
        -40176 * t100 * t115
        + 324 * t100 * t25
        - 8154 * t100 * t35
        - 11664 * t100 * t5
        + 2592 * t103 * t21
        + 40176 * t106 * t108
        + 1944 * t106 * t12
        + 15984 * t120 * t115
        + 216 * t120 * t25
        + 2916 * t120 * t35
        + 5184 * t120 * t5
    )
    t144 = t3 * numberOfSamples
    t166 = (
        4 * t144 * powerSums["43"]
        + t144 * powerSums["61"]
        - 2600 * t15 * powerSums["43"]
        - 650 * t15 * powerSums["61"]
        - 68 * t3 * powerSums["43"]
        - 17 * t3 * powerSums["61"]
        + 5631 * t38 * powerSums["25"]
        + 7508 * t38 * powerSums["43"]
        + 1877 * t38 * powerSums["61"]
        + 560 * t9 * powerSums["43"]
        + 140 * t9 * powerSums["61"]
    )
    t191 = (
        2808 * numberOfSamples * powerSums["25"]
        + 3744 * numberOfSamples * powerSums["43"]
        + 936 * numberOfSamples * powerSums["61"]
        - 9612 * t1 * powerSums["25"]
        - 12816 * t1 * powerSums["43"]
        - 3204 * t1 * powerSums["61"]
        - 10851 * t2 * powerSums["25"]
        - 14468 * t2 * powerSums["43"]
        - 3617 * t2 * powerSums["61"]
        + 13602 * t8 * powerSums["25"]
        + 18136 * t8 * powerSums["43"]
        + 4534 * t8 * powerSums["61"]
    )
    t201 = t10 ** 2
    t208 = t9 * powerSums["23"]
    t210 = t9 * powerSums["41"]
    t212 = t15 * powerSums["23"]
    t214 = t15 * powerSums["41"]
    t216 = t38 * powerSums["23"]
    t218 = t38 * powerSums["41"]
    t220 = t2 * powerSums["23"]
    t222 = t2 * powerSums["41"]
    t225 = t8 * powerSums["11"]
    t226 = t225 * powerSums["12"]
    t228 = t8 * powerSums["23"]
    t230 = t8 * powerSums["41"]
    t232 = t1 * powerSums["11"]
    t233 = t232 * powerSums["12"]
    t235 = t1 * powerSums["23"]
    t237 = t1 * powerSums["41"]
    t239 = numberOfSamples * powerSums["11"]
    t240 = t239 * powerSums["12"]
    t242 = numberOfSamples * powerSums["23"]
    t244 = numberOfSamples * powerSums["41"]
    t246 = powerSums["11"] * powerSums["12"]
    t248 = powerSums["11"] * powerSums["30"]
    t250 = (
        14982 * t226
        + 6388 * t228
        + 1494 * t230
        - 34893 * t233
        - 6696 * t235
        - 1944 * t237
        + 43362 * t240
        + 2664 * t242
        + 864 * t244
        - 22680 * t246
        - 7560 * t248
    )
    t252 = t15 * powerSums["11"]
    t253 = t252 * powerSums["12"]
    t255 = t15 * powerSums["20"]
    t256 = t255 * powerSums["21"]
    t258 = t38 * powerSums["11"]
    t259 = t258 * powerSums["12"]
    t261 = t38 * powerSums["20"]
    t262 = t261 * powerSums["21"]
    t264 = t2 * powerSums["11"]
    t265 = t264 * powerSums["12"]
    t267 = t2 * powerSums["20"]
    t268 = t267 * powerSums["21"]
    t270 = t8 * powerSums["20"]
    t271 = t270 * powerSums["21"]
    t273 = t232 * powerSums["30"]
    t275 = t1 * powerSums["20"]
    t276 = t275 * powerSums["21"]
    t278 = t239 * powerSums["30"]
    t280 = numberOfSamples * powerSums["20"]
    t281 = t280 * powerSums["21"]
    t283 = (
        78 * t253
        + 27 * t256
        + 57 * t259
        - 57 * t262
        - 3057 * t265
        - 117 * t268
        + 600 * t271
        - 11667 * t273
        - 774 * t276
        + 14526 * t278
        + 324 * t281
    )
    t284 = t225 * powerSums["30"]
    t286 = t264 * powerSums["30"]
    t288 = t258 * powerSums["30"]
    t290 = t252 * powerSums["30"]
    t292 = 5 * t38
    t301 = t9 * powerSums["20"]
    t304 = t9 * powerSums["11"]
    t307 = numberOfSamples * powerSums["03"]
    t317 = (-2 + numberOfSamples) ** 2
    t318 = 11 * t8
    t328 = 7 * t15
    t336 = 5 * t2
    t353 = 4 * t38
    t360 = t10 * powerSums["21"]
    t363 = (
        4866 * t284
        - 891 * t286
        - 21 * t288
        + 30 * t290
        - 9
        * (t292 - 45 * t2 + 89 * t8 + 293 * t1 - 1182 * numberOfSamples + 840)
        * t16
        * powerSums["11"]
        - 3 * t301 * powerSums["21"]
        - 3 * t304 * powerSums["30"]
        - t307
        * (
            -3 * (3 * t38 - 37 * t2 + 195 * t8 - 539 * t1 + 834 * numberOfSamples - 696) * t10
            + t317 * (t2 - t318 + 67 * t1 - 165 * numberOfSamples + 108) * powerSums["20"]
        )
        - 9 * t304 * powerSums["12"]
        + t120
        * (
            (-t328 + 75 * t38 - 289 * t2 + 321 * t8 + 512 * t1 - 1332 * numberOfSamples + 1440)
            * powerSums["13"]
            + 18
            * (t38 - t336 - 42 * t8 + 375 * t1 - 947 * numberOfSamples + 738)
            * powerSums["11"]
            * powerSums["20"]
            + (
                -t328
                + 87 * t38
                - 445 * t2
                + 1173 * t8
                - 1780 * t1
                + 1548 * numberOfSamples
                + 144
            )
            * powerSums["31"]
        )
        + 6
        * numberOfSamples
        * (t353 - 49 * t2 + 248 * t8 - 659 * t1 + 1032 * numberOfSamples - 936)
        * t360
    )
    t368 = powerSums["10"] * powerSums["12"]
    t371 = t1 * powerSums["03"]
    t374 = t8 * powerSums["03"]
    t377 = t2 * powerSums["03"]
    t380 = (
        -1950 * t15 * powerSums["25"]
        + 420 * t9 * powerSums["25"]
        - 51 * t3 * powerSums["25"]
        + 3 * t144 * powerSums["25"]
        - 22680 * t201 * powerSums["21"]
        - 7560 * powerSums["03"] * t201
        - 3
        * powerSums["02"]
        * (
            t3 * powerSums["23"]
            + t3 * powerSums["41"]
            - 2 * t208
            - 12 * t210
            - 80 * t212
            + 52 * t214
            + 742 * t216
            - 42 * t218
            - 3017 * t220
            - 413 * t222
            + t250
            + t283
            + t363
        )
        + 46602 * t307 * t368
        - 49365 * t371 * t368
        + 30462 * t374 * t368
        - 11808 * t377 * t368
    )
    t381 = t9 * powerSums["03"]
    t384 = t15 * powerSums["03"]
    t385 = t10 * powerSums["20"]
    t390 = t38 * powerSums["03"]
    t399 = powerSums["10"] * powerSums["30"]
    t410 = (
        7884 * t307 * t385
        + 13374 * t307 * t399
        + 24 * t381 * t368
        - 387 * t384 * t368
        + 2832 * t390 * t368
        - 11178 * t371 * t385
        - 12231 * t371 * t399
        + 7812 * t374 * t385
        + 7618 * t374 * t399
        - 2826 * t377 * t385
        - 36 * t384 * t385
        + 504 * t390 * t385
    )
    t415 = t9 * powerSums["04"]
    t416 = powerSums["10"] * powerSums["11"]
    t419 = powerSums["11"] * powerSums["22"]
    t422 = powerSums["12"] * powerSums["21"]
    t425 = powerSums["13"] * powerSums["20"]
    t428 = powerSums["20"] * powerSums["31"]
    t431 = powerSums["11"] * t34
    t434 = powerSums["11"] * powerSums["20"]
    t439 = t15 * powerSums["04"]
    t444 = (
        36 * t24 * t419
        + 54 * t24 * t422
        + 21 * t24 * t425
        + 21 * t24 * t428
        - 891 * t33 * t422
        - 27 * t33 * t431
        + 10 * t381 * t399
        - 159 * t384 * t399
        + 171 * t39 * t434
        + 6 * t415 * t416
        - 36 * t439 * t416
    )
    t455 = t38 * powerSums["04"]
    t470 = (
        -261 * t33 * t425
        - 225 * t33 * t428
        - 3522 * t377 * t399
        + 1030 * t390 * t399
        + 6 * t455 * t416
        + 6858 * t52 * t422
        + 1335 * t52 * t425
        + 867 * t52 * t428
        + 81 * t52 * t431
        + 1674 * t67 * t431
        - 1575 * t55 * t434
        + 3663 * t70 * t434
    )
    t472 = t2 * powerSums["04"]
    t481 = powerSums["20"] * powerSums["21"]
    t488 = t8 * powerSums["04"]
    t499 = (
        306 * t472 * t416
        - 552 * t488 * t416
        - 30510 * t67 * t422
        + 83970 * t84 * t422
        - 3519 * t67 * t425
        + 5340 * t84 * t425
        - 963 * t67 * t428
        - 1536 * t84 * t428
        - 12015 * t84 * t431
        + 6111 * t87 * t434
        + 18909 * t81 * t481
        - 27756 * t96 * t481
    )
    t504 = t1 * powerSums["04"]
    t517 = numberOfSamples * powerSums["04"]
    t526 = (
        -142911 * t100 * t422
        - 4644 * t100 * t425
        + 3996 * t100 * t428
        + 28269 * t100 * t431
        - 31050 * t103 * t434
        + 21060 * t106 * t481
        + 138510 * t120 * t422
        - 432 * t120 * t425
        - 4320 * t120 * t428
        - 21222 * t120 * t431
        + 162 * t504 * t416
        + 108 * t517 * t416
    )
    t530 = t3 * powerSums["03"]
    t532 = t3 * powerSums["11"]
    t539 = powerSums["20"] * powerSums["30"]
    t554 = (
        72 * t252 * t21
        + 9 * t304 * t21
        + 36 * t252 * t539
        - 99 * t30 * t481
        + 3 * t304 * t539
        + 60 * t304 * powerSums["14"]
        + 3 * t381 * t34
        - 45 * t384 * t34
        + 16 * t381 * powerSums["40"]
        - t530 * powerSums["40"]
        - 6 * t532 * powerSums["14"]
    )
    t579 = (
        -2277 * t258 * t21
        + 16641 * t264 * t21
        - 150 * t252 * powerSums["14"]
        - 879 * t258 * t539
        - 714 * t258 * powerSums["14"]
        + 5931 * t264 * t539
        - 1053 * t377 * t34
        + 291 * t390 * t34
        - 112 * t384 * powerSums["40"]
        + 478 * t390 * powerSums["40"]
        + 1287 * t49 * t481
        - 6921 * t64 * t481
    )
    t585 = powerSums["11"] * powerSums["40"]
    t588 = t8 * t107
    t603 = t1 * t107
    t606 = (
        162 * t100 * t585
        - 59238 * t225 * t21
        - 20130 * t225 * t539
        - 15546 * t225 * powerSums["14"]
        + 5772 * t264 * powerSums["14"]
        + 9786 * t374 * powerSums["22"]
        + 3106 * t374 * powerSums["40"]
        - 1471 * t377 * powerSums["40"]
        - 552 * t84 * t585
        + 6606 * t588 * powerSums["21"]
        - 6588 * t603 * powerSums["21"]
    )
    t619 = numberOfSamples * t107
    t632 = (
        108 * t120 * t585
        + 118503 * t232 * t21
        - 135270 * t239 * t21
        + 39393 * t232 * t539
        + 18576 * t232 * powerSums["14"]
        - 44874 * t239 * t539
        - 7992 * t239 * powerSums["14"]
        + 2808 * t307 * powerSums["22"]
        + 1584 * t307 * powerSums["40"]
        - 8316 * t371 * powerSums["22"]
        - 3600 * t371 * powerSums["40"]
        + 2592 * t619 * powerSums["21"]
    )
    t637 = t9 * powerSums["05"]
    t640 = t9 * t107
    t645 = t15 * powerSums["05"]
    t648 = t15 * t107
    t653 = t38 * powerSums["05"]
    t656 = t38 * t107
    t661 = t2 * powerSums["05"]
    t664 = (
        3 * t637 * t10
        - 36 * t645 * t10
        + 210 * t653 * t10
        - 612 * t661 * t10
        + 156 * t252 * powerSums["32"]
        - 2520 * t258 * powerSums["32"]
        + 36 * t304 * powerSums["32"]
        - 6 * t532 * powerSums["32"]
        + 18 * t640 * powerSums["21"]
        - 216 * t648 * powerSums["21"]
        + 1152 * t656 * powerSums["21"]
    )
    t665 = t2 * t107
    t672 = t8 * powerSums["05"]
    t679 = t1 * powerSums["05"]
    t686 = numberOfSamples * powerSums["05"]
    t693 = (
        759 * t672 * t10
        - 216 * t679 * t10
        - 108 * t686 * t10
        - 25164 * t225 * powerSums["32"]
        + 26352 * t232 * powerSums["32"]
        - 10368 * t239 * powerSums["32"]
        + 11514 * t264 * powerSums["32"]
        + 612 * t661 * powerSums["20"]
        - 3564 * t665 * powerSums["21"]
        - 759 * t672 * powerSums["20"]
        + 216 * t679 * powerSums["20"]
        + 108 * t686 * powerSums["20"]
    )
    t719 = (
        -1944 * t307 * t107
        + 4680 * t371 * t107
        - 3990 * t374 * t107
        + 1404 * t377 * t107
        - 5961 * t377 * powerSums["22"]
        + 48 * t381 * powerSums["22"]
        - 408 * t384 * powerSums["22"]
        + 2046 * t390 * powerSums["22"]
        - 3 * t530 * powerSums["22"]
        - 3 * t637 * powerSums["20"]
        + 36 * t645 * powerSums["20"]
        - 210 * t653 * powerSums["20"]
    )
    t726 = powerSums["02"] ** 2
    t727 = numberOfSamples * t726
    t761 = (
        -120 * t390 * t107
        - 36 * t384 * t107
        + 6 * t381 * t107
        - 9
        * t727
        * (
            3
            * (t38 - 3 * t2 - 62 * t8 + 445 * t1 - 1047 * numberOfSamples + 786)
            * powerSums["10"]
            * powerSums["11"]
            + (-t292 + 56 * t2 - 238 * t8 + 478 * t1 - 453 * numberOfSamples + 162)
            * powerSums["21"]
        )
        - 108 * t517 * powerSums["21"]
        - 162 * t504 * powerSums["21"]
        + 552 * t488 * powerSums["21"]
        - 306 * t472 * powerSums["21"]
        - 6 * t455 * powerSums["21"]
        + 36 * t439 * powerSums["21"]
        - 6 * t415 * powerSums["21"]
        + 10350 * t307 * t201
    )
    t766 = t3 * powerSums["13"]
    t769 = t3 * powerSums["20"]
    t772 = t3 * powerSums["21"]
    t775 = t3 * powerSums["30"]
    t778 = t9 * t34
    t783 = t9 * powerSums["13"]
    t788 = t9 * powerSums["21"]
    t791 = t9 * powerSums["30"]
    t794 = t15 * t34
    t797 = (
        33 * t390 * t201
        + 12 * t301 * powerSums["23"]
        - 5 * t766 * powerSums["30"]
        - 3 * t769 * powerSums["23"]
        - 3 * t772 * powerSums["40"]
        - 3 * t775 * powerSums["31"]
        + 9 * t778 * powerSums["21"]
        + 91 * t783 * powerSums["30"]
        + 54 * t788 * powerSums["40"]
        + 47 * t791 * powerSums["31"]
        - 126 * t794 * powerSums["21"]
    )
    t800 = t15 * powerSums["13"]
    t805 = t15 * powerSums["21"]
    t808 = t15 * powerSums["30"]
    t811 = t38 * t34
    t816 = t38 * powerSums["13"]
    t821 = t38 * powerSums["21"]
    t824 = t38 * powerSums["30"]
    t827 = t2 * t34
    t830 = (
        2505 * t374 * t201
        - 447 * t377 * t201
        + 168 * t255 * powerSums["23"]
        - 1806 * t261 * powerSums["23"]
        - 716 * t800 * powerSums["30"]
        - 444 * t805 * powerSums["40"]
        - 348 * t808 * powerSums["31"]
        + 729 * t811 * powerSums["21"]
        + 3154 * t816 * powerSums["30"]
        + 2052 * t821 * powerSums["40"]
        + 1574 * t824 * powerSums["31"]
        - 2169 * t827 * powerSums["21"]
    )
    t834 = t2 * powerSums["13"]
    t839 = t2 * powerSums["21"]
    t842 = t2 * powerSums["30"]
    t847 = t8 * t34
    t850 = t8 * powerSums["13"]
    t855 = t8 * powerSums["21"]
    t858 = t8 * powerSums["30"]
    t861 = (
        -7041 * t371 * t201
        + 7827 * t267 * powerSums["23"]
        - 17646 * t270 * powerSums["23"]
        + 2190 * t374 * t34
        - 8693 * t834 * powerSums["30"]
        - 5655 * t839 * powerSums["40"]
        - 4815 * t842 * powerSums["31"]
        + 3420 * t847 * powerSums["21"]
        + 15259 * t850 * powerSums["30"]
        + 9234 * t855 * powerSums["40"]
        + 9683 * t858 * powerSums["31"]
    )
    t864 = t1 * powerSums["13"]
    t869 = t1 * powerSums["21"]
    t872 = t1 * powerSums["30"]
    t877 = numberOfSamples * powerSums["13"]
    t882 = numberOfSamples * powerSums["21"]
    t885 = numberOfSamples * powerSums["30"]
    t888 = powerSums["03"] * powerSums["10"]
    t893 = (
        19656 * t275 * powerSums["23"]
        - 8208 * t280 * powerSums["23"]
        + 972 * t307 * t34
        - 2358 * t371 * t34
        - 15498 * t864 * powerSums["30"]
        - 8154 * t869 * powerSums["40"]
        - 11106 * t872 * powerSums["31"]
        + 7848 * t877 * powerSums["30"]
        + 2916 * t882 * powerSums["40"]
        + 6408 * t885 * powerSums["31"]
        - 22680 * t888 * powerSums["12"]
        - 7560 * t888 * powerSums["30"]
    )
    t898 = t3 * powerSums["12"]
    t911 = t9 * powerSums["12"]
    t920 = (
        -51 * t17 * powerSums["31"]
        + 244 * t24 * powerSums["33"]
        - 12 * t304 * powerSums["50"]
        + 645 * t39 * powerSums["31"]
        - 16 * t4 * powerSums["33"]
        - 9 * t772 * powerSums["22"]
        + 144 * t788 * powerSums["22"]
        - 15 * t898 * powerSums["13"]
        - 9 * t898 * powerSums["31"]
        + 261 * t911 * powerSums["13"]
        + 153 * t911 * powerSums["31"]
    )
    t925 = t15 * powerSums["12"]
    t938 = t38 * powerSums["12"]
    t947 = (
        162 * t252 * powerSums["50"]
        - 966 * t258 * powerSums["50"]
        - 1804 * t33 * powerSums["33"]
        + 7324 * t52 * powerSums["33"]
        - 3567 * t55 * powerSums["31"]
        + 10083 * t70 * powerSums["31"]
        - 1152 * t805 * powerSums["22"]
        + 5526 * t821 * powerSums["22"]
        - 2004 * t925 * powerSums["13"]
        - 1188 * t925 * powerSums["31"]
        + 8838 * t938 * powerSums["13"]
        + 5346 * t938 * powerSums["31"]
    )
    t955 = t2 * powerSums["12"]
    t968 = t8 * powerSums["12"]
    t975 = (
        -6582 * t225 * powerSums["50"]
        + 3294 * t264 * powerSums["50"]
        - 17476 * t67 * powerSums["33"]
        - 6036 * t67 * powerSums["51"]
        - 16335 * t839 * powerSums["22"]
        + 26128 * t84 * powerSums["33"]
        + 8532 * t84 * powerSums["51"]
        + 28890 * t855 * powerSums["22"]
        - 25143 * t955 * powerSums["13"]
        - 15381 * t955 * powerSums["31"]
        + 46413 * t968 * powerSums["13"]
        + 28413 * t968 * powerSums["31"]
    )
    t982 = t1 * powerSums["12"]
    t995 = numberOfSamples * powerSums["12"]
    t1002 = (
        -23328 * t100 * powerSums["33"]
        - 7128 * t100 * powerSums["51"]
        + 8928 * t120 * powerSums["33"]
        + 2592 * t120 * powerSums["51"]
        + 6912 * t232 * powerSums["50"]
        - 2808 * t239 * powerSums["50"]
        - 27432 * t869 * powerSums["22"]
        + 10368 * t882 * powerSums["22"]
        - 49302 * t982 * powerSums["13"]
        - 30510 * t982 * powerSums["31"]
        + 25272 * t995 * powerSums["13"]
        + 17496 * t995 * powerSums["31"]
    )
    t1018 = powerSums["21"] * powerSums["30"]
    t1027 = t107 * powerSums["11"]
    t1030 = (
        24 * t24 * t1018
        - 36 * t33 * t1027
        + 30 * t11 * powerSums["23"]
        + 21 * t11 * powerSums["41"]
        - 39 * t17 * powerSums["13"]
        + 6 * t24 * t585
        - 117 * t30 * t246
        - 45 * t30 * t248
        - 6 * t4 * powerSums["15"]
        - 6 * t4 * powerSums["51"]
        - 3 * t769 * powerSums["41"]
    )
    t1031 = t38 * t201
    t1034 = t201 * powerSums["10"]
    t1035 = t2 * t1034
    t1058 = (
        -387 * t33 * t1018
        + 108 * t1031 * powerSums["21"]
        - 135 * t1035 * powerSums["11"]
        + 90 * t24 * powerSums["15"]
        + 90 * t24 * powerSums["51"]
        + 1359 * t49 * t246
        - 396 * t30 * powerSums["23"]
        - 288 * t30 * powerSums["41"]
        + 33 * t301 * powerSums["41"]
        - 324 * t33 * t419
        - 36 * t33 * t585
        + 489 * t39 * powerSums["13"]
    )
    t1064 = t2 * t201
    t1067 = t8 * t1034
    t1084 = (
        360 * t52 * t1027
        - 1404 * t1064 * powerSums["21"]
        + 1458 * t1067 * powerSums["11"]
        + 459 * t49 * t248
        - 120 * t255 * powerSums["41"]
        - 654 * t33 * powerSums["15"]
        - 654 * t33 * powerSums["51"]
        + 1008 * t52 * t419
        + 2460 * t49 * powerSums["23"]
        + 1830 * t49 * powerSums["41"]
        + 6 * t52 * t585
    )
    t1095 = t8 * t201
    t1098 = t1 * t1034
    t1111 = (
        2832 * t52 * t1018
        - 1368 * t67 * t1027
        + 7668 * t1095 * powerSums["21"]
        - 5589 * t1098 * powerSums["11"]
        - 6741 * t64 * t246
        - 1845 * t64 * t248
        - 84 * t261 * powerSums["41"]
        + 2610 * t52 * powerSums["15"]
        + 2610 * t52 * powerSums["51"]
        - 2715 * t55 * powerSums["13"]
        - 7920 * t64 * powerSums["23"]
        - 6084 * t64 * powerSums["41"]
    )
    t1128 = t1 * t201
    t1131 = numberOfSamples * t1034
    t1138 = (
        -11808 * t67 * t1018
        + 2448 * t84 * t1027
        - 21276 * t1128 * powerSums["21"]
        + 7506 * t1131 * powerSums["11"]
        + 17649 * t81 * t246
        + 3861 * t81 * t248
        + 1851 * t267 * powerSums["41"]
        - 216 * t67 * t419
        + 306 * t67 * t585
        - 6036 * t67 * powerSums["15"]
        + 7791 * t70 * powerSums["13"]
    )
    t1157 = numberOfSamples * t201
    t1164 = (
        -2052 * t100 * t1027
        + 30462 * t84 * t1018
        + 31104 * t1157 * powerSums["21"]
        - 25758 * t96 * t246
        - 5238 * t96 * t248
        - 5241 * t270 * powerSums["41"]
        - 5580 * t84 * t419
        + 12630 * t81 * powerSums["23"]
        + 10353 * t81 * powerSums["41"]
        + 8532 * t84 * powerSums["15"]
        - 11574 * t87 * powerSums["13"]
        - 14454 * t87 * powerSums["31"]
    )
    t1174 = t1 * t34
    t1187 = t16 * powerSums["11"]
    t1192 = (
        -49365 * t100 * t1018
        + 10908 * t100 * t419
        - 7128 * t100 * powerSums["15"]
        + 648 * t120 * t1027
        + 8208 * t103 * powerSums["13"]
        + 9504 * t103 * powerSums["31"]
        + 20088 * t106 * t246
        + 4968 * t106 * t248
        - 2673 * t1174 * powerSums["21"]
        + 22680 * t1187 * powerSums["20"]
        - 8964 * t96 * powerSums["23"]
        - 8316 * t96 * powerSums["41"]
    )
    t1203 = numberOfSamples * t34
    t1219 = (
        -22680 * powerSums["10"] * powerSums["21"] * powerSums["30"]
        + 46602 * t120 * t1018
        + 2160 * t106 * powerSums["23"]
        + 2484 * t106 * powerSums["41"]
        - 5832 * t120 * t419
        + 2592 * t120 * powerSums["15"]
        + 810 * t1203 * powerSums["21"]
        + 68040 * t246 * powerSums["20"]
        + 6048 * t275 * powerSums["41"]
        - 2484 * t280 * powerSums["41"]
        - 68040 * t368 * powerSums["21"]
        + 22680 * t434 * powerSums["30"]
    )
    t1229 = powerSums["13"] ** 2
    t1242 = t107 ** 2
    t1257 = (
        2592 * numberOfSamples * t1242
        - 6264 * t1 * t1242
        - t144 * t1229
        - 878 * t15 * t1229
        - 18065 * t2 * t1229
        + 7 * t3 * t1229
        + 5191 * t38 * t1229
        + 38102 * t8 * t1229
        + 52 * t9 * t1229
        - 2340 * t2 * t1242
        + 5580 * t8 * t1242
    )
    t1259 = t3 * t1
    t1281 = (
        -36 * t15 * t1242
        + 468 * t38 * t1242
        + 2 * t1259 * powerSums["44"]
        + t1259 * powerSums["62"]
        - 34 * t144 * powerSums["44"]
        - 17 * t144 * powerSums["62"]
        + 3754 * t15 * powerSums["44"]
        + 280 * t3 * powerSums["44"]
        + 140 * t3 * powerSums["62"]
        - 1300 * t9 * powerSums["44"]
        - 650 * t9 * powerSums["62"]
    )
    t1296 = powerSums["31"] ** 2
    t1305 = (
        18720 * numberOfSamples * t1296
        - 41688 * t1 * t1296
        + 1872 * t1 * powerSums["44"]
        + 936 * t1 * powerSums["62"]
        + 1877 * t15 * powerSums["62"]
        + 9068 * t2 * powerSums["44"]
        + 4534 * t2 * powerSums["62"]
        - 7234 * t38 * powerSums["44"]
        - 3617 * t38 * powerSums["62"]
        - 6408 * t8 * powerSums["44"]
        - 3204 * t8 * powerSums["62"]
    )
    t1325 = powerSums["22"] ** 2
    t1351 = (
        33696 * numberOfSamples * t1325
        - 82296 * t1 * t1325
        + 936 * t1 * powerSums["26"]
        - 18065 * t2 * t1296
        + 5191 * t38 * t1296
        + 38102 * t8 * t1296
        - 38106 * t2 * t1325
        + 77256 * t8 * t1325
        + 4534 * t2 * powerSums["26"]
        - 3617 * t38 * powerSums["26"]
        - 3204 * t8 * powerSums["26"]
    )
    t1359 = powerSums["01"] ** 2
    t1360 = t1359 ** 2
    t1361 = 2 * t8
    t1364 = t1361 - 21 * t1 + 79 * numberOfSamples - 105
    t1369 = 207 * numberOfSamples
    t1370 = 5 * t8 - 54 * t1 + t1369 - 278
    t1374 = numberOfSamples * powerSums["02"]
    t1382 = -2 * t2 + t318 - 16 * t1 + numberOfSamples + 6
    t1387 = 13 * t2
    t1394 = 11 * t38
    t1410 = 6 * t1
    t1412 = t8 - t1410 + 11 * numberOfSamples - 6
    t1413 = -4 + numberOfSamples
    t1416 = 5 * numberOfSamples
    t1417 = t1 - t1416 + 16
    t1420 = 7 * numberOfSamples
    t1440 = 7 * t2
    t1462 = (
        -7560 * t360
        + powerSums["03"]
        * (
            (17 * t38 - 191 * t2 + 937 * t8 - 2449 * t1 + 3486 * numberOfSamples - 2520) * t10
            - numberOfSamples
            * (t38 - t1440 + 35 * t8 - 125 * t1 + 204 * numberOfSamples - 108)
            * powerSums["20"]
        )
        + 288 * t244
        - 900 * t237
        + 1066 * t230
        - 608 * t222
        + 180 * t218
        - 28 * t214
        + 2 * t210
        + 216 * t242
        - 1044 * t235
        + 1572 * t228
        - 1016 * t220
        + 320 * t216
        - 52 * t212
    )
    t1474 = -1 + numberOfSamples
    t1484 = (-3 + numberOfSamples) ** 2
    t1486 = 6 * numberOfSamples
    t1505 = (
        4 * t208
        - 10638 * t120 * t434
        + 2637 * t100 * t434
        + 801 * t84 * t434
        - 405 * t67 * t434
        + 45 * t52 * t434
        + 3
        * t1474
        * powerSums["02"]
        * (
            (19 * t2 - 156 * t8 + 251 * t1 + 930 * numberOfSamples - 2520)
            * powerSums["10"]
            * powerSums["11"]
            - t1484 * numberOfSamples * (t1 + t1486 - 16) * powerSums["21"]
        )
        + 6696 * t240
        - 8514 * t233
        + 5679 * t226
        - 2043 * t265
        + 369 * t259
        - 27 * t253
        - 531 * t64 * powerSums["21"]
        + 45 * t49 * powerSums["21"]
        + 2736 * t120 * powerSums["31"]
    )
    t1531 = (
        -17 * t33 * powerSums["13"]
        - 13 * t33 * powerSums["31"]
        + 215 * t52 * powerSums["13"]
        + 163 * t52 * powerSums["31"]
        - 102 * t55 * powerSums["11"]
        - 1189 * t67 * powerSums["13"]
        - 905 * t67 * powerSums["31"]
        + 1044 * t70 * powerSums["11"]
        + 2709 * t81 * powerSums["21"]
        - 6 * t256
        + 60 * t262
        - 228 * t268
        - 551 * t286
        + 97 * t288
        - 7 * t290
    )
    t1556 = t416 * powerSums["20"]
    t1558 = (
        -4818 * t100 * powerSums["13"]
        - 3858 * t100 * powerSums["31"]
        + 5076 * t103 * powerSums["11"]
        + 10422 * t106 * powerSums["21"]
        + 3168 * t120 * powerSums["13"]
        + 3361 * t84 * powerSums["13"]
        + 2597 * t84 * powerSums["31"]
        - 3858 * t87 * powerSums["11"]
        - 7245 * t96 * powerSums["21"]
        + 7560 * t1556
        + 408 * t271
        - 2490 * t273
        - 342 * t276
        + 2088 * t278
        + 108 * t281
        + 1583 * t284
    )
    t1565 = t201 * t10
    t1583 = powerSums["12"] ** 2
    t1585 = powerSums["12"] * powerSums["30"]
    t1587 = (
        -5688 * numberOfSamples * t1565
        + 1512 * t1 * t1565
        - 20916 * t103 * powerSums["30"]
        - 144 * t8 * t1565
        + 45360 * t16 * powerSums["12"]
        + 15120 * t16 * powerSums["30"]
        - 135 * t49 * t34
        + 1431 * t64 * t34
        + 14694 * t87 * powerSums["30"]
        + 7560 * t1565
        + 68040 * t1583
        + 45360 * t1585
    )
    t1620 = 2 * t38
    t1631 = 4 * numberOfSamples
    t1652 = 35 * t1
    t1663 = powerSums["30"] ** 2
    t1664 = numberOfSamples * t1663
    t1666 = t1 * t1663
    t1668 = t8 * t1663
    t1670 = t2 * t1663
    t1672 = (
        -5622 * t70 * powerSums["30"]
        + 1146 * t55 * powerSums["30"]
        - 102 * t39 * powerSums["30"]
        + 216 * t120 * powerSums["50"]
        + 432 * t100 * powerSums["50"]
        - 6
        * t1374
        * (
            -9 * t1370 * t201
            - 3 * (t353 - 25 * t2 - 94 * t8 + 1093 * t1 - 2850 * numberOfSamples + 2232) * t107
            + 9
            * (t336 - 56 * t8 + 219 * t1 - 300 * numberOfSamples + 12)
            * t10
            * powerSums["20"]
            - 3
            * powerSums["10"]
            * (
                (t1394 - 143 * t2 + 769 * t8 - 2101 * t1 + 3084 * numberOfSamples - 2340)
                * powerSums["12"]
                + 2
                * (t1620 - 28 * t2 + 157 * t8 - 434 * t1 + 621 * numberOfSamples - 438)
                * powerSums["30"]
            )
            + t1412
            * (
                -3 * (t1 - t1631 - 3) * t34
                + 3 * (t8 - t1 - t1486 + 18) * powerSums["22"]
                + (t8 - t1410 + t1416 + 24) * powerSums["40"]
            )
        )
        + 9
        * t727
        * (
            -3 * (t336 - 53 * t8 + 201 * t1 - 267 * numberOfSamples - 6) * t10
            + (t38 - t1440 + 20 * t8 - t1652 + 39 * numberOfSamples - 18) * powerSums["20"]
        )
        + 297 * t455 * powerSums["20"]
        - 27 * t439 * powerSums["20"]
        - 14670 * t1664
        + 14103 * t1666
        - 7666 * t1668
        + 2604 * t1670
    )
    t1674 = t3 * powerSums["60"]
    t1675 = t9 * t1583
    t1677 = powerSums["21"] ** 2
    t1678 = t9 * t1677
    t1680 = t9 * t1663
    t1682 = t34 * powerSums["20"]
    t1684 = 9 * t15 * t1682
    t1685 = t15 * t1583
    t1687 = t15 * t1677
    t1689 = t15 * t1663
    t1692 = 117 * t38 * t1682
    t1693 = t38 * t1583
    t1695 = t38 * t1677
    t1697 = t38 * t1663
    t1699 = (
        t1674
        - 54 * t1675
        - 18 * t1678
        - 4 * t1680
        + t1684
        + 927 * t1685
        + 234 * t1687
        + 69 * t1689
        - t1692
        - 6894 * t1693
        - 1386 * t1695
        - 556 * t1697
    )
    t1701 = 612 * t2 * t1682
    t1702 = t2 * t1583
    t1704 = t2 * t1677
    t1707 = 1557 * t8 * t1682
    t1708 = t8 * t1583
    t1710 = t8 * t1677
    t1713 = 1863 * t1 * t1682
    t1714 = t1 * t1583
    t1716 = t1 * t1677
    t1719 = 810 * numberOfSamples * t1682
    t1720 = numberOfSamples * t1583
    t1722 = numberOfSamples * t1677
    t1725 = 936 * numberOfSamples * powerSums["60"]
    t1726 = (
        t1701
        + 28746 * t1702
        + 4518 * t1704
        - t1707
        - 75726 * t1708
        - 8100 * t1710
        + t1713
        + 129951 * t1714
        + 7344 * t1716
        - t1719
        - 132030 * t1720
        - 2592 * t1722
        - t1725
    )
    t1729 = t3 * powerSums["24"]
    t1731 = t3 * powerSums["42"]
    t1743 = t9 * powerSums["24"]
    t1745 = t9 * powerSums["42"]
    t1748 = 16 * t9 * powerSums["60"]
    t1753 = (
        -42 * t24 * powerSums["14"]
        - 6 * t24 * powerSums["50"]
        + 144 * t30 * powerSums["22"]
        + 27 * t30 * powerSums["40"]
        + 576 * t33 * powerSums["14"]
        + 72 * t33 * powerSums["50"]
        - 270 * t39 * powerSums["12"]
        + 15 * t1729
        + 12 * t1731
        - 228 * t1743
        - 180 * t1745
        - t1748
    )
    t1760 = t15 * powerSums["24"]
    t1762 = t15 * powerSums["42"]
    t1765 = 124 * t15 * powerSums["60"]
    t1776 = t38 * powerSums["24"]
    t1778 = t38 * powerSums["42"]
    t1780 = (
        -1656 * t49 * powerSums["22"]
        - 297 * t49 * powerSums["40"]
        - 3660 * t52 * powerSums["14"]
        - 420 * t52 * powerSums["50"]
        + 3186 * t55 * powerSums["12"]
        + 8712 * t64 * powerSums["22"]
        + 1539 * t64 * powerSums["40"]
        - 16254 * t70 * powerSums["12"]
        + 1680 * t1760
        + 1308 * t1762
        + t1765
        - 6798 * t1776
        - 5220 * t1778
    )
    t1783 = 526 * t38 * powerSums["60"]
    t1794 = t2 * powerSums["24"]
    t1796 = t2 * powerSums["42"]
    t1799 = 1351 * t2 * powerSums["60"]
    t1808 = (
        12168 * t67 * powerSums["14"]
        + 1224 * t67 * powerSums["50"]
        - 23832 * t81 * powerSums["22"]
        - 4239 * t81 * powerSums["40"]
        - 20706 * t84 * powerSums["14"]
        - 1518 * t84 * powerSums["50"]
        + 43470 * t87 * powerSums["12"]
        + 30888 * t96 * powerSums["22"]
        + 5562 * t96 * powerSums["40"]
        - t1783
        + 16125 * t1794
        + 12072 * t1796
        + t1799
    )
    t1811 = t8 * powerSums["24"]
    t1813 = t8 * powerSums["42"]
    t1816 = 2266 * t8 * powerSums["60"]
    t1823 = t1 * powerSums["24"]
    t1825 = t1 * powerSums["42"]
    t1828 = 2268 * t1 * powerSums["60"]
    t1831 = numberOfSamples * powerSums["24"]
    t1833 = numberOfSamples * powerSums["42"]
    t1835 = (
        16632 * t100 * powerSums["14"]
        - 62532 * t103 * powerSums["12"]
        - 14256 * t106 * powerSums["22"]
        - 2592 * t106 * powerSums["40"]
        - 4968 * t120 * powerSums["14"]
        - 23862 * t1811
        - 17064 * t1813
        - t1816
        + 21060 * t1823
        + 14256 * t1825
        + t1828
        - 7992 * t1831
        - 5184 * t1833
    )
    t1839 = t381 * powerSums["21"]
    t1843 = t304 * powerSums["13"]
    t1845 = t304 * powerSums["31"]
    t1847 = t911 * powerSums["30"]
    t1849 = t301 * powerSums["22"]
    t1852 = 6 * t301 * powerSums["40"]
    t1863 = (
        27 * t439 * t10
        + 234 * t33 * t18
        + 144 * t33 * t21
        - 60 * t24 * powerSums["32"]
        + 54 * t33 * t539
        + 90 * t384 * t416
        - 6 * t1839
        - 54 * t1843
        - 30 * t1845
        - 30 * t1847
        - 18 * t1849
        - t1852
    )
    t1864 = t648 * powerSums["20"]
    t1870 = t384 * powerSums["21"]
    t1874 = t252 * powerSums["13"]
    t1876 = t252 * powerSums["31"]
    t1878 = t925 * powerSums["30"]
    t1880 = t255 * powerSums["22"]
    t1883 = 81 * t255 * powerSums["40"]
    t1890 = (
        -297 * t455 * t10
        + 297 * t1064 * powerSums["20"]
        - 378 * t49 * t107
        - 2718 * t52 * t18
        + 792 * t33 * powerSums["32"]
        - 918 * t390 * t416
        + 36 * t1864
        + 54 * t1870
        + 666 * t1874
        + 306 * t1876
        + 540 * t1878
        + 144 * t1880
        + t1883
    )
    t1896 = t656 * powerSums["20"]
    t1902 = t390 * powerSums["21"]
    t1906 = t258 * powerSums["13"]
    t1908 = t258 * powerSums["31"]
    t1910 = t938 * powerSums["30"]
    t1912 = t261 * powerSums["22"]
    t1915 = 483 * t261 * powerSums["40"]
    t1918 = (
        3780 * t64 * t107
        - 3078 * t1095 * powerSums["20"]
        - 1764 * t52 * t21
        + 3690 * t377 * t416
        - 666 * t52 * t539
        - 4920 * t52 * powerSums["32"]
        + 18 * t1896
        - 438 * t1902
        - 3294 * t1906
        - 1110 * t1908
        - 4278 * t1910
        - 324 * t1912
        - t1915
    )
    t1927 = t665 * powerSums["20"]
    t1935 = t377 * powerSums["21"]
    t1941 = t264 * powerSums["13"]
    t1943 = t264 * powerSums["31"]
    t1945 = (
        1539 * t472 * t10
        - 13662 * t81 * t107
        + 11475 * t1128 * powerSums["20"]
        + 13482 * t67 * t18
        + 8928 * t67 * t21
        - 5427 * t81 * t34
        - 1539 * t472 * powerSums["20"]
        + 3510 * t67 * t539
        + 15840 * t67 * powerSums["32"]
        - 3924 * t1927
        + 2682 * t1935
        + 7974 * t1941
        + 990 * t1943
    )
    t1948 = t955 * powerSums["30"]
    t1950 = t267 * powerSums["22"]
    t1953 = 1485 * t267 * powerSums["40"]
    t1964 = t588 * powerSums["20"]
    t1972 = (
        -4239 * t488 * t10
        + 18036 * t96 * t107
        - 15174 * t1157 * powerSums["20"]
        - 35298 * t84 * t18
        - 23724 * t84 * t21
        + 7209 * t96 * t34
        - 7722 * t374 * t416
        - 9702 * t84 * t539
        + 19098 * t1948
        - 144 * t1950
        + t1953
        + 24606 * t1964
    )
    t1973 = t374 * powerSums["21"]
    t1979 = t225 * powerSums["13"]
    t1981 = t225 * powerSums["31"]
    t1983 = t968 * powerSums["30"]
    t1985 = t270 * powerSums["22"]
    t1988 = 2319 * t270 * powerSums["40"]
    t1999 = (
        5562 * t504 * t10
        + 51516 * t100 * t18
        + 37152 * t100 * t21
        + 15012 * t100 * t539
        + 10476 * t371 * t416
        + 4239 * t488 * powerSums["20"]
        - 25260 * t84 * powerSums["32"]
        - 8196 * t1973
        - 10044 * t1979
        + 2436 * t1981
        - 52824 * t1983
        + 1422 * t1985
        - t1988
    )
    t2001 = t603 * powerSums["20"]
    t2007 = t371 * powerSums["21"]
    t2013 = t232 * powerSums["13"]
    t2015 = t232 * powerSums["31"]
    t2017 = t982 * powerSums["30"]
    t2019 = t275 * powerSums["22"]
    t2022 = 1674 * t275 * powerSums["40"]
    t2027 = (
        -2592 * t517 * t10
        + 17928 * t100 * powerSums["32"]
        - 1296 * t106 * t107
        + 162 * t106 * t34
        - 9936 * t307 * t416
        - 5562 * t504 * powerSums["20"]
        - 56376 * t2001
        + 11088 * t2007
        + 6480 * t2013
        - 5184 * t2015
        + 91386 * t2017
        - 1728 * t2019
        + t2022
    )
    t2034 = t619 * powerSums["20"]
    t2037 = 5184 * t307 * powerSums["21"]
    t2042 = t239 * powerSums["13"]
    t2044 = t239 * powerSums["31"]
    t2046 = t995 * powerSums["30"]
    t2048 = t280 * powerSums["22"]
    t2051 = 432 * t280 * powerSums["40"]
    t2053 = (
        -40176 * t120 * t18
        - 33696 * t120 * t21
        - 12528 * t120 * t539
        - 4320 * t120 * powerSums["32"]
        + 2592 * t517 * powerSums["20"]
        + 7560 * t1663
        + 42120 * t2034
        - t2037
        + 2592 * t2042
        + 6912 * t2044
        - 90612 * t2046
        + 648 * t2048
        - t2051
    )
    t2062 = numberOfSamples * t1412
    t2069 = powerSums["03"] ** 2
    t2099 = 9 * t1
    t2145 = t10 - powerSums["20"]
    t2158 = (
        9 * numberOfSamples * t1382 * t201 * powerSums["20"]
        + 6
        * t2062
        * t16
        * (
            6 * (t1 - t1631 + 3) * powerSums["12"]
            + (t1 - numberOfSamples + 18) * powerSums["30"]
        )
        - 3
        * numberOfSamples
        * (2 * t15 - 27 * t38 + 161 * t2 - 495 * t8 + 773 * t1 - 558 * numberOfSamples + 144)
        * powerSums["04"]
        * t2145
        + t3 * t1663
        - 1152 * t1664
        + 3564 * t1666
        - 4420 * t1668
        + 2851 * t1670
        - 1030 * t1697
        + 208 * t1689
        - 22 * t1680
        + 136080 * t108
        + t1719
    )
    t2168 = (
        54 * t1678
        - t1684
        - 918 * t1687
        + t1692
        + 6534 * t1695
        - t1701
        - 24570 * t1704
        + t1707
        + 50436 * t1710
        - t1713
        - 15876 * t1714
        - 52272 * t1716
        + 7776 * t1720
        + 20736 * t1722
    )
    t2182 = (
        9 * t3 * t1583
        - t1674
        - 144 * t1675
        + 882 * t1685
        - 2232 * t1693
        - 24 * t1731
        + 420 * t1743
        + 372 * t1745
        - 2952 * t1760
        - 2580 * t1762
        + 12030 * t1776
        + 10452 * t1778
        + t1783
    )
    t2193 = (
        297 * t1702
        - 30393 * t1794
        - 26340 * t1796
        - t1799
        + 9288 * t1708
        + 46950 * t1811
        + 40152 * t1813
        + t1816
        - 40068 * t1823
        - 33264 * t1825
        - t1828
        + 14040 * t1831
        + 11232 * t1833
        + t1725
    )
    t2244 = (
        -t317
        * numberOfSamples
        * (t38 - t1440 + 8 * t8 + 88 * t1 - 297 * numberOfSamples + 207)
        * powerSums["14"]
        - t9 * powerSums["50"]
        + 36 * numberOfSamples * powerSums["50"]
        + 72 * t1 * powerSums["50"]
        - 253 * t8 * powerSums["50"]
        + 204 * t2 * powerSums["50"]
        - 70 * t38 * powerSums["50"]
        + 12 * t15 * powerSums["50"]
        + 22680 * t18
        - t3 * powerSums["32"]
        - 2736 * numberOfSamples * powerSums["32"]
        + 6552 * t1 * powerSums["32"]
        - 5882 * t8 * powerSums["32"]
        + 2609 * t2 * powerSums["32"]
        - 602 * t38 * powerSums["32"]
        + 56 * t15 * powerSums["32"]
        + 4 * t9 * powerSums["32"]
        + (
            t9
            + 12 * t15
            - 293 * t38
            + 1977 * t2
            - 6710 * t8
            + 13131 * t1
            - 14958 * numberOfSamples
            + 7560
        )
        * powerSums["03"]
        * powerSums["11"]
        - 45090 * t239 * powerSums["21"]
    )
    t2284 = (
        -19746 * t225 * powerSums["21"]
        + 39501 * t232 * powerSums["21"]
        + 24 * t252 * powerSums["21"]
        - 15 * t255 * powerSums["30"]
        - 759 * t258 * powerSums["21"]
        + 115 * t261 * powerSums["30"]
        + 5547 * t264 * powerSums["21"]
        - 477 * t267 * powerSums["30"]
        + 1036 * t270 * powerSums["30"]
        - 1092 * t275 * powerSums["30"]
        + 432 * t280 * powerSums["30"]
        + t301 * powerSums["30"]
        + 3 * t304 * powerSums["21"]
        + 3 * t911 * powerSums["20"]
        - 27 * t925 * powerSums["20"]
        + 57 * t938 * powerSums["20"]
        + 117 * t955 * powerSums["20"]
        - 600 * t968 * powerSums["20"]
        + 774 * t982 * powerSums["20"]
        - 324 * t995 * powerSums["20"]
    )
    t2296 = (
        -27 * t1729
        - 6 * powerSums["10"] * (t2244 + t2284)
        - t1765
        + t1748
        + 96 * t1843
        + 14256 * t2044
        - 15228 * t2015
        + 11682 * t1981
        - 10170 * t1943
        + 4662 * t1908
        - 954 * t1876
        + 72 * t1845
        + t2037
    )
    t2330 = (
        -12096 * t2007
        + 10692 * t1973
        - 4950 * t1935
        + 1386 * t1902
        - 234 * t1870
        + 18 * t1839
        + 18576 * t2042
        - 26892 * t2013
        + 24162 * t1979
        - 17154 * t1941
        + 6846 * t1906
        - 1314 * t1874
        - 9
        * t106
        * (
            -2 * (t1620 + t2 - 218 * t8 + 1367 * t1 - 3132 * numberOfSamples + 2340) * t107
            + t1412
            * (
                (-t1 + numberOfSamples - 3) * t34
                + 2 * (t8 - 2 * t1 - t1416 + 6) * powerSums["22"]
                + 3 * t1417 * powerSums["40"]
            )
        )
        + 9 * t2062 * t1565
    )
    t2344 = (
        6 * t898 * powerSums["30"]
        - 114 * t1847
        + 36 * t1849
        + t1852
        - 90 * t1864
        + 942 * t1878
        - 324 * t1880
        - t1883
        + 522 * t1896
        - 4290 * t1910
        + 252 * t1912
        + t1915
        + 5112 * t1927
        + 11364 * t1948
    )
    t2355 = (
        5724 * t1950
        - t1953
        - 55494 * t1964
        - 17196 * t1983
        - 22320 * t1985
        + t1988
        + 188298 * t2001
        + 13608 * t2017
        + 30888 * t2019
        - t2022
        - 267948 * t2034
        - 4320 * t2046
        - 14256 * t2048
        + t2051
    )
    t2360 = (
        -2268 * t15 * t1325
        + 288 * t9 * t1325
        - 18 * t3 * t1325
        + 3
        * t1360
        * (
            -48 * t1364 * t201
            + 18 * numberOfSamples * t1370 * t385
            + 3
            * t1374
            * (
                3 * (t318 - 114 * t1 + 425 * numberOfSamples - 562) * t10
                + t1382 * powerSums["20"]
            )
            - 2
            * powerSums["10"]
            * (
                36
                * (t38 - t1387 + 71 * t8 - 197 * t1 + 288 * numberOfSamples - 210)
                * powerSums["12"]
                + (t1394 - 149 * t2 + 835 * t8 - 2347 * t1 + 3450 * numberOfSamples - 2520)
                * powerSums["30"]
            )
            + numberOfSamples
            * (
                (-39 * t2 + 444 * t8 - 1761 * t1 + 2400 * numberOfSamples + 36) * t107
                + t1412
                * (
                    -6 * t1413 * t34
                    + 9 * t1417 * powerSums["22"]
                    + 2 * (t1 - t1420 + 18) * powerSums["40"]
                )
            )
        )
        - 6 * t1359 * powerSums["01"] * (t1462 + t1505 + t1531 + t1558)
        + t1359
        * (
            t1587
            + t1672
            + t1699
            + t1726
            + t1753
            + t1780
            + t1808
            + t1835
            + t1863
            + t1890
            + t1918
            + t1945
            + t1972
            + t1999
            + t2027
            + t2053
        )
        - 9 * t1360 * t1359 * (8 * t1364 * t10 - t2062 * powerSums["20"])
        + 68040 * t10 * t1677
        + 7560 * t2069 * t10
        - 3
        * t726
        * (
            6 * numberOfSamples * (t2 - 10 * t8 + t1652 - 50 * numberOfSamples + 24) * t201
            + 3
            * (
                t328
                - 55 * t38
                - 160 * t2
                + 2809 * t8
                - 10179 * t1
                + 14778 * numberOfSamples
                - 7560
            )
            * t107
            - 6
            * numberOfSamples
            * (t38 - 10 * t2 + 32 * t8 - 32 * t1 - 9 * numberOfSamples + 18)
            * t385
            + 2
            * t2062
            * powerSums["10"]
            * (
                3 * (t8 - 8 * t1 + 22 * numberOfSamples - 15) * powerSums["12"]
                + (t8 - t2099 + 32 * numberOfSamples - 54) * powerSums["30"]
            )
            - t2062
            * (
                -6 * (t1 - t1420 + 9) * t34
                + 3 * (t1361 - 5 * t1 - 45 * numberOfSamples + 144) * powerSums["22"]
                + 2 * (t8 - t2099 + 26 * numberOfSamples - 30) * powerSums["40"]
            )
        )
        + powerSums["02"] * (t2158 + t2168 + t2182 + t2193 + t2296 + t2330 + t2344 + t2355)
    )
    t2361 = powerSums["11"] * powerSums["31"]
    t2370 = powerSums["11"] * powerSums["23"]
    t2385 = (
        -52704 * t100 * t2370
        + 2592 * t106 * t2361
        - 54 * t11 * t2361
        + 20736 * t120 * t2370
        + 666 * t30 * t2361
        - 3294 * t49 * t2361
        + 7974 * t64 * t2361
        - 10044 * t81 * t2361
        + 6480 * t96 * t2361
        - 23028 * t67 * t2370
        + 50328 * t84 * t2370
    )
    t2390 = powerSums["12"] * powerSums["22"]
    t2395 = t107 * powerSums["30"]
    t2402 = powerSums["11"] * powerSums["41"]
    t2407 = powerSums["22"] * powerSums["30"]
    t2435 = (
        -450 * t49 * t108
        - 1692 * t64 * t108
        - 438 * t49 * t1585
        - 312 * t33 * t2370
        + 5040 * t52 * t2370
        + 2304 * t33 * t2390
        + 72 * t33 * t2395
        + 240 * t52 * t2395
        + 300 * t33 * t2402
        + 1428 * t52 * t2402
        + 816 * t33 * t2407
    )
    t2459 = (
        19674 * t81 * t108
        - 51300 * t96 * t108
        + 2682 * t64 * t1585
        - 8196 * t81 * t1585
        - 2808 * t67 * t2395
        + 7980 * t84 * t2395
        - 11544 * t67 * t2402
        + 31092 * t84 * t2402
        - 4092 * t52 * t2407
        + 11922 * t67 * t2407
        - 19572 * t84 * t2407
    )
    t2482 = (
        54864 * t100 * t2390
        - 9360 * t100 * t2395
        - 37152 * t100 * t2402
        + 16632 * t100 * t2407
        - 5184 * t106 * t1585
        - 20736 * t120 * t2390
        + 3888 * t120 * t2395
        + 15984 * t120 * t2402
        - 5616 * t120 * t2407
        + 11088 * t96 * t1585
        - 136080 * t416 * t481
    )
    t2487 = powerSums["10"] * powerSums["13"]
    t2496 = powerSums["13"] * powerSums["21"]
    t2499 = powerSums["21"] * powerSums["31"]
    t2506 = powerSums["11"] * powerSums["13"]
    t2509 = t107 * powerSums["12"]
    t2535 = (
        42 * t384 * t1187
        - 306 * t24 * t2496
        - 522 * t24 * t2499
        - 60 * t381 * t246
        - 94 * t381 * t2487
        + 306 * t30 * t2506
        + 432 * t33 * t2509
        + 540 * t384 * t360
        - 108 * t415 * t368
        - 114 * t381 * t481
        + 24 * t637 * t416
    )
    t2580 = (
        3306 * t377 * t1187
        - 11052 * t52 * t2390
        + 3088 * t390 * t248
        - 3148 * t390 * t2487
        - 10692 * t52 * t2496
        - 17676 * t52 * t2499
        + 990 * t64 * t2506
        - 2304 * t52 * t2509
        + 19098 * t377 * t360
        - 4104 * t455 * t368
        + 1932 * t653 * t416
    )
    t2612 = powerSums["10"] * powerSums["31"]
    t2627 = (
        14940 * t371 * t1187
        - 57780 * t84 * t2390
        + 30904 * t374 * t248
        - 19366 * t374 * t2487
        - 56826 * t84 * t2496
        - 92826 * t84 * t2499
        + 2436 * t81 * t2506
        - 13212 * t84 * t2509
        - 30518 * t374 * t2612
        - 18468 * t488 * t368
        + 13164 * t672 * t416
    )
    t2651 = (
        61020 * t100 * t2496
        + 13176 * t100 * t2509
        - 40176 * t103 * t18
        - 12528 * t307 * t1187
        - 33048 * t371 * t248
        + 22212 * t371 * t2487
        - 5184 * t96 * t2506
        + 30996 * t371 * t2612
        + 91386 * t371 * t360
        + 16308 * t504 * t368
        - 13824 * t679 * t416
    )
    t2674 = (
        98604 * t100 * t2499
        + 6912 * t106 * t2506
        - 34992 * t120 * t2496
        - 50544 * t120 * t2499
        - 5184 * t120 * t2509
        + 13248 * t307 * t248
        - 12816 * t307 * t2487
        - 15696 * t307 * t2612
        - 90612 * t307 * t360
        - 5832 * t517 * t368
        + 5616 * t686 * t416
    )
    t2698 = t3 * powerSums["04"]
    t2728 = (
        2 * t2698 * t399
        - 312 * t455 * t385
        + 468 * t472 * t385
        + 318 * t488 * t385
        - 1404 * t504 * t385
        + 864 * t517 * t385
        - 32 * t415 * t399
        + 224 * t439 * t399
        - 956 * t455 * t399
        + 2942 * t472 * t399
        - 6212 * t488 * t399
    )
    t2752 = (
        -954 * t252 * t425
        - 1314 * t252 * t428
        + 4662 * t258 * t425
        + 6846 * t258 * t428
        - 10170 * t264 * t425
        + 72 * t304 * t425
        + 11364 * t377 * t481
        + 942 * t384 * t481
        - 6 * t415 * t385
        + 72 * t439 * t385
        - 4290 * t390 * t481
    )
    t2775 = (
        11682 * t225 * t425
        + 24162 * t225 * t428
        - 15228 * t232 * t425
        - 26892 * t232 * t428
        + 14256 * t239 * t425
        + 18576 * t239 * t428
        - 17154 * t264 * t428
        - 4320 * t307 * t481
        + 13608 * t371 * t481
        - 17196 * t374 * t481
        - 12096 * t982 * t539
    )
    t2821 = (
        24 * t252 * t1018
        + 2244 * t258 * t1018
        + 2244 * t390 * t246
        + 1980 * t252 * t422
        - 6152 * t252 * powerSums["33"]
        - 2580 * t252 * powerSums["51"]
        - 6840 * t258 * t422
        + 872 * t304 * powerSums["33"]
        + 96 * t415 * powerSums["22"]
        - 816 * t439 * powerSums["22"]
        - 234 * t925 * t539
    )
    t2845 = (
        33240 * t225 * t1018
        - 13206 * t264 * t1018
        + 216 * t225 * t422
        + 33240 * t374 * t246
        - 13206 * t377 * t246
        + 10452 * t258 * powerSums["51"]
        + 10314 * t264 * t422
        - 26340 * t264 * powerSums["51"]
        + 4524 * t455 * powerSums["22"]
        + 1386 * t938 * t539
        - 4950 * t955 * t539
    )
    t2868 = (
        -38664 * t232 * t1018
        + 16416 * t239 * t1018
        + 40152 * t225 * powerSums["51"]
        - 15768 * t232 * t422
        - 33264 * t232 * powerSums["51"]
        + 10368 * t239 * t422
        + 11232 * t239 * powerSums["51"]
        + 16416 * t307 * t246
        - 38664 * t371 * t246
        + 10692 * t968 * t539
        + 5184 * t995 * t539
    )
    t2874 = t3 * powerSums["06"]
    t2883 = t9 * powerSums["06"]
    t2899 = t15 * powerSums["06"]
    t2912 = t38 * powerSums["06"]
    t2919 = (
        124 * t2899 * t10
        - 526 * t2912 * t10
        + 426 * t472 * t201
        - 1146 * t488 * t201
        - 124 * t2899 * powerSums["20"]
        + 372 * t304 * powerSums["51"]
        - 90 * t439 * t34
        + 546 * t455 * t34
        - 212 * t439 * powerSums["40"]
        + 800 * t455 * powerSums["40"]
        + 324 * t645 * powerSums["21"]
    )
    t2927 = t2 * powerSums["06"]
    t2947 = t8 * powerSums["06"]
    t2956 = t1 * powerSums["06"]
    t2965 = numberOfSamples * powerSums["06"]
    t2968 = (
        2266 * t2947 * powerSums["20"]
        - 2268 * t2956 * powerSums["20"]
        + 936 * t2965 * powerSums["20"]
        - 2916 * t504 * t34
        + 1080 * t517 * t34
        + 3920 * t488 * powerSums["40"]
        - 4320 * t504 * powerSums["40"]
        + 1872 * t517 * powerSums["40"]
        - 13164 * t672 * powerSums["21"]
        + 13824 * t679 * powerSums["21"]
        - 5616 * t686 * powerSums["21"]
    )
    t2998 = t144 * powerSums["03"]
    t3006 = t1474 ** 2
    t3025 = (
        1804 * t384 * powerSums["41"]
        - 352 * t381 * powerSums["41"]
        + 40 * t530 * powerSums["41"]
        - 2 * t2998 * powerSums["41"]
        - 18
        * numberOfSamples
        * t1360
        * powerSums["01"]
        * (
            -3 * t1370 * powerSums["10"] * powerSums["11"]
            + t3006 * (t1 - t1416 + 6) * powerSums["21"]
        )
        - 8208 * t517 * t107
        + 18360 * t504 * t107
        - 14280 * t488 * t107
        + 4608 * t472 * t107
        - 420 * t455 * t107
        - 72 * t439 * t107
    )
    t3029 = t9 * t2069
    t3036 = t15 * t2069
    t3041 = t38 * t2069
    t3046 = t2 * t2069
    t3053 = (
        t2069 * t3 * powerSums["20"]
        - 4 * t3029 * t10
        + 69 * t3036 * t10
        - 556 * t3041 * t10
        + 2604 * t3046 * t10
        + 12 * t415 * t107
        - 22 * t3029 * powerSums["20"]
        + 208 * t3036 * powerSums["20"]
        - 1030 * t3041 * powerSums["20"]
        + 2851 * t3046 * powerSums["20"]
        - 10006 * t390 * powerSums["23"]
    )
    t3054 = t8 * t2069
    t3061 = t1 * t2069
    t3068 = numberOfSamples * t2069
    t3080 = (
        45360 * t10 * powerSums["03"] * powerSums["21"]
        - 7666 * t3054 * t10
        + 14103 * t3061 * t10
        - 14670 * t3068 * t10
        - 4420 * t3054 * powerSums["20"]
        + 3564 * t3061 * powerSums["20"]
        - 1152 * t3068 * powerSums["20"]
        - 18720 * t307 * powerSums["23"]
        + 48096 * t371 * powerSums["23"]
        - 49208 * t374 * powerSums["23"]
        + 27748 * t377 * powerSums["23"]
    )
    t3086 = t144 * powerSums["21"]
    t3091 = t434 * powerSums["21"]
    t3132 = (
        -816 * t15 * powerSums["22"] * powerSums["40"]
        + 96 * t9 * powerSums["22"] * powerSums["40"]
        - 126 * t390 * t1556
        + 468 * t33 * t3091
        + 342 * t52 * t3091
        - 400 * t381 * powerSums["23"]
        + 2452 * t384 * powerSums["23"]
        - 828 * t648 * powerSums["22"]
        - 72 * t648 * powerSums["40"]
        - 420 * t656 * powerSums["40"]
        - 1080 * t788 * powerSums["41"]
    )
    t3158 = (
        -15918 * t2 * powerSums["22"] * powerSums["40"]
        + 4524 * t38 * powerSums["22"] * powerSums["40"]
        - 70002 * t371 * t1556
        + 29196 * t374 * t1556
        - 5346 * t377 * t1556
        + 25112 * t258 * powerSums["33"]
        - 63488 * t264 * powerSums["33"]
        - 18342 * t67 * t3091
        + 89892 * t84 * t3091
        - 14280 * t588 * powerSums["40"]
        + 4608 * t665 * powerSums["40"]
    )
    t3184 = (
        14040 * numberOfSamples * powerSums["22"] * powerSums["40"]
        - 34668 * t1 * powerSums["22"] * powerSums["40"]
        + 32748 * t8 * powerSums["22"] * powerSums["40"]
        - 209358 * t100 * t3091
        + 260172 * t120 * t3091
        + 87156 * t307 * t1556
        + 98432 * t225 * powerSums["33"]
        - 84672 * t232 * powerSums["33"]
        + 29952 * t239 * powerSums["33"]
        + 75168 * t869 * powerSums["41"]
        - 28080 * t882 * powerSums["41"]
    )
    t3187 = t144 * powerSums["12"]
    t3232 = (
        -2580 * t252 * powerSums["15"]
        + 10452 * t258 * powerSums["15"]
        + 3276 * t656 * powerSums["22"]
        - 24 * t783 * powerSums["31"]
        - 908 * t800 * powerSums["31"]
        + 6708 * t805 * powerSums["23"]
        + 7182 * t816 * powerSums["31"]
        + 5736 * t925 * powerSums["14"]
        + 6708 * t925 * powerSums["32"]
        - 20358 * t938 * powerSums["14"]
        - 26154 * t938 * powerSums["32"]
    )
    t3256 = (
        40152 * t225 * powerSums["15"]
        - 26340 * t264 * powerSums["15"]
        - 6804 * t588 * powerSums["22"]
        - 3708 * t665 * powerSums["22"]
        - 26154 * t821 * powerSums["23"]
        - 27770 * t834 * powerSums["31"]
        + 70068 * t839 * powerSums["23"]
        + 50304 * t955 * powerSums["14"]
        + 70068 * t955 * powerSums["32"]
        - 81804 * t968 * powerSums["14"]
        - 121296 * t968 * powerSums["32"]
    )
    t3279 = (
        -33264 * t232 * powerSums["15"]
        + 11232 * t239 * powerSums["15"]
        + 18360 * t603 * powerSums["22"]
        - 10368 * t619 * powerSums["22"]
        + 60524 * t850 * powerSums["31"]
        - 66096 * t864 * powerSums["31"]
        + 29952 * t877 * powerSums["31"]
        + 75168 * t982 * powerSums["14"]
        + 116640 * t982 * powerSums["32"]
        - 28080 * t995 * powerSums["14"]
        - 44928 * t995 * powerSums["32"]
    )
    t3328 = (
        882 * t255 * t1677
        - 2232 * t261 * t1677
        - 918 * t1685 * powerSums["20"]
        + 6534 * t1693 * powerSums["20"]
        - 24570 * t1702 * powerSums["20"]
        + 495 * t656 * t34
        + 1440 * t665 * t34
        + 5736 * t805 * powerSums["41"]
        + 2452 * t808 * powerSums["32"]
        - 20358 * t821 * powerSums["41"]
        - 10006 * t824 * powerSums["32"]
    )
    t3373 = (
        16686 * t1174 * powerSums["22"]
        - 7776 * t1203 * powerSums["22"]
        - 15876 * t275 * t1677
        + 7776 * t280 * t1677
        - 52272 * t1714 * powerSums["20"]
        + 20736 * t1720 * powerSums["20"]
        - 133002 * t619 * t34
        + 116640 * t869 * powerSums["23"]
        + 48096 * t872 * powerSums["32"]
        - 44928 * t882 * powerSums["23"]
        - 18720 * t885 * powerSums["32"]
    )
    t3376 = t144 * powerSums["10"]
    t3382 = t9 * t16
    t3425 = (
        1804 * t15 * powerSums["14"] * powerSums["30"]
        - 352 * t9 * powerSums["14"] * powerSums["30"]
        - 2580 * t255 * powerSums["24"]
        + 372 * t301 * powerSums["24"]
        + 1300 * t33 * powerSums["16"]
        - 1080 * t39 * powerSums["14"]
        - 3754 * t52 * powerSums["16"]
        + 3648 * t55 * powerSums["14"]
        - 153 * t794 * powerSums["22"]
        + 63 * t811 * powerSums["22"]
        + 324 * t925 * powerSums["50"]
    )
    t3451 = (
        14572 * t2 * powerSums["14"] * powerSums["30"]
        - 6142 * t38 * powerSums["14"] * powerSums["30"]
        + 10452 * t261 * powerSums["24"]
        - 26340 * t267 * powerSums["24"]
        - 6396 * t70 * powerSums["14"]
        + 3123 * t827 * powerSums["22"]
        - 11961 * t847 * powerSums["22"]
        + 5400 * t87 * powerSums["14"]
        - 1932 * t938 * powerSums["50"]
        + 6588 * t955 * powerSums["50"]
        - 13164 * t968 * powerSums["50"]
    )
    t3477 = (
        -7488 * numberOfSamples * powerSums["14"] * powerSums["30"]
        + 20448 * t1 * powerSums["14"] * powerSums["30"]
        - 22880 * t8 * powerSums["14"] * powerSums["30"]
        - 1728 * t103 * powerSums["14"]
        + 40152 * t270 * powerSums["24"]
        - 33264 * t275 * powerSums["24"]
        + 11232 * t280 * powerSums["24"]
        + 18360 * t603 * powerSums["40"]
        - 8208 * t619 * powerSums["40"]
        + 13824 * t982 * powerSums["50"]
        - 5616 * t995 * powerSums["50"]
    )
    t3485 = t3 * t10
    t3527 = (
        -297 * t1031 * powerSums["22"]
        + 126 * t1035 * powerSums["12"]
        + 1332 * t1064 * t107
        - 180 * t11 * powerSums["24"]
        - 228 * t11 * powerSums["42"]
        + 234 * t30 * t1583
        + 927 * t30 * t1677
        + 312 * t17 * powerSums["32"]
        + 136 * t4 * powerSums["34"]
        + 102 * t4 * powerSums["52"]
        - 27 * t769 * powerSums["42"]
    )
    t3551 = (
        1539 * t1064 * powerSums["22"]
        - 306 * t1067 * powerSums["12"]
        - 5283 * t1095 * t107
        - 1386 * t49 * t1583
        - 6894 * t49 * t1677
        - 1120 * t24 * powerSums["34"]
        - 840 * t24 * powerSums["52"]
        + 1308 * t30 * powerSums["24"]
        + 1680 * t30 * powerSums["42"]
        + 420 * t301 * powerSums["42"]
        - 1920 * t39 * powerSums["32"]
    )
    t3574 = (
        7200 * t1128 * t107
        - 4239 * t1095 * powerSums["22"]
        + 306 * t1098 * powerSums["12"]
        + 4518 * t64 * t1583
        + 28746 * t64 * t1677
        - 2952 * t255 * powerSums["42"]
        + 5200 * t33 * powerSums["34"]
        + 3900 * t33 * powerSums["52"]
        - 5220 * t49 * powerSums["24"]
        - 6798 * t49 * powerSums["42"]
        + 6096 * t55 * powerSums["32"]
    )
    t3620 = (
        -2266 * t2947 * t10
        + 108 * t1157 * t107
        + 7344 * t96 * t1583
        + 129951 * t96 * t1677
        - 30393 * t267 * powerSums["42"]
        + 7234 * t67 * powerSums["16"]
        + 28936 * t67 * powerSums["34"]
        + 21702 * t67 * powerSums["52"]
        - 17064 * t81 * powerSums["24"]
        - 23862 * t81 * powerSums["42"]
        + 6264 * t87 * powerSums["32"]
    )
    t3644 = (
        2268 * t2956 * t10
        - 1296 * t103 * powerSums["32"]
        - 2592 * t106 * t1583
        - 132030 * t106 * t1677
        - 2592 * t1157 * powerSums["22"]
        + 46950 * t270 * powerSums["42"]
        - 9068 * t84 * powerSums["16"]
        - 36272 * t84 * powerSums["34"]
        - 27204 * t84 * powerSums["52"]
        + 14256 * t96 * powerSums["24"]
        + 21060 * t96 * powerSums["42"]
    )
    t3667 = (
        -936 * t2965 * t10
        + 6408 * t100 * powerSums["16"]
        + 25632 * t100 * powerSums["34"]
        + 19224 * t100 * powerSums["52"]
        - 5184 * t106 * powerSums["24"]
        - 7992 * t106 * powerSums["42"]
        - 1872 * t120 * powerSums["16"]
        - 7488 * t120 * powerSums["34"]
        - 5616 * t120 * powerSums["52"]
        - 40068 * t275 * powerSums["42"]
        + 14040 * t280 * powerSums["42"]
    )
    uglyD1O4E = (
        (
            t3256
            + t3232
            + t3184
            + t3158
            + t3132
            + t3080
            + t3053
            + t3025
            + t2968
            + t2919
            + t2868
            + t2845
            + t2821
            + t2775
            + t2752
            - 8100 * t81 * t1583
            - 75726 * t81 * t1677
            - 18 * t11 * t1583
            - 54 * t11 * t1677
            - 117 * t1031 * t107
            + 297 * t267 * t1677
            - 25281 * t588 * t34
            + 9288 * t270 * t1677
            + 91611 * t603 * t34
            + 9 * t769 * t1677
            - 144 * t301 * t1677
            - 63 * t648 * t34
            - 18 * t381 * t1556
            - 54 * t24 * t3091
            + 180 * t384 * t1556
            - 1746 * t472 * t34
            + 1351 * t2927 * t10
            + 1440 * t504 * t201
            + 3120 * t488 * t34
            - 648 * t517 * t201
            + t2874 * t10
            + 6 * t439 * t201
            + 6 * t415 * t34
            - 16 * t2883 * t10
            - 78 * t455 * t201
            + 2 * t530 * t248
            + 18 * t532 * t422
            + 6 * t532 * t1018
            - 8 * t381 * t248
            - 288 * t304 * t422
            + 96 * t304 * t428
            - 60 * t304 * t1018
            + 18 * t911 * t539
            + t2728
            + t2674
            + t2651
            + t2627
            + t2580
            + t2535
            + t2482
            + t2459
            + t2435
            + t2385
            + t2360
            + 10 * t530 * t2612
            - 45360 * t888 * t434
            + 6 * t2698 * t368
            - 3168 * t517 * t399
            + 7200 * t504 * t399
            + 17386 * t377 * t2612
            - 6308 * t390 * t2612
            + 1432 * t384 * t2612
            - 182 * t381 * t2612
            + 7128 * t67 * t2509
            - 9498 * t374 * t1187
            + 9630 * t377 * t2487
            - 13918 * t377 * t248
            + 11310 * t472 * t368
            - 6588 * t661 * t416
            + 32670 * t67 * t2390
            + 30762 * t67 * t2496
            + 50286 * t67 * t2499
            - 52824 * t374 * t360
            - 1110 * t49 * t2506
            - 582 * t390 * t1187
            + 696 * t384 * t2487
            + 24 * t384 * t246
            - 268 * t384 * t248
            + 888 * t439 * t368
            - 324 * t645 * t416
            + 2376 * t33 * t2496
            + 4008 * t33 * t2499
            - 4278 * t390 * t360
            + 6 * t4 * t2407
            - 30 * t381 * t360
            - 30 * t11 * t2506
            - 36 * t24 * t2509
            + 6 * t530 * t2487
            + 6 * t530 * t246
            + 6 * t530 * t481
            + 12 * t4 * t2402
            + 18 * t4 * t2496
            + 30 * t4 * t2499
            - 6 * t11 * t1585
            - 12 * t24 * t2395
            + 72 * t30 * t108
            - 72 * t24 * t2370
            - 120 * t24 * t2402
            - 288 * t24 * t2390
            - 96 * t24 * t2407
            + 54 * t30 * t1585
            + 12 * t4 * t2370
            + 18 * t4 * t2390
            + 9
            * numberOfSamples
            * (t38 - t1387 + 68 * t8 - 173 * t1 + t1369 - 90)
            * t726
            * powerSums["02"]
            * t2145
            - t144 * t1296
            + 7 * t3 * t1296
            + 52 * t9 * t1296
            - 878 * t15 * t1296
            + 11448 * t38 * t1325
            + t1351
            + t1305
            + t1281
            + t1257
            + 68040 * t107 * t34
            - 41688 * t1 * t1229
            - 6 * t3 * powerSums["22"] * powerSums["40"]
            - 2 * t144 * powerSums["13"] * powerSums["31"]
            - 2 * t144 * powerSums["30"] * powerSums["32"]
            - 2 * t144 * powerSums["14"] * powerSums["30"]
            + 40 * t3 * powerSums["14"] * powerSums["30"]
            + 27 * t15 * t201 * powerSums["22"]
            - 18 * t38 * t1034 * powerSums["12"]
            - 90 * t33 * t35
            - 2214 * t39 * t18
            - 18 * t39 * t21
            - 2
            * powerSums["01"]
            * (
                t1219
                + t1192
                + t1164
                + t1138
                + t1111
                + t1084
                + t1058
                + t1030
                + t1002
                + t975
                + t947
                + t920
                + t893
                + t861
                + t830
                + t797
                + t761
                + t719
                + t693
                + t664
                + t632
                + t606
                + t579
                + t554
                + t526
                + t499
                + t470
                + t444
                + t191
                + t166
                + t380
                + t410
            )
            + 6 * t4 * t5
            - 18 * t11 * t12
            + 162 * t17 * t18
            + 18 * t17 * t21
            + t1259 * powerSums["26"]
            - 17 * t144 * powerSums["26"]
            + 140 * t3 * powerSums["26"]
            + 12 * t24 * t25
            - 650 * t9 * powerSums["26"]
            + 1877 * t15 * powerSums["26"]
            + 18720 * numberOfSamples * t1229
            - 72 * t24 * t5
            + 126 * t30 * t12
            - 8 * t3376 * powerSums["34"]
            - 6 * t3376 * powerSums["52"]
            + 12 * t3485 * powerSums["24"]
            + 15 * t3485 * powerSums["42"]
            - 24 * t3382 * powerSums["32"]
            - 2 * t3376 * powerSums["16"]
            - 12 * t3382 * powerSums["14"]
            + 34 * t4 * powerSums["16"]
            - 24 * t769 * powerSums["24"]
            + 18 * t778 * powerSums["22"]
            + 168 * t17 * powerSums["14"]
            - 280 * t24 * powerSums["16"]
            - 24 * t911 * powerSums["50"]
            - 400 * t791 * powerSums["32"]
            + 50304 * t839 * powerSums["41"]
            + 27748 * t842 * powerSums["32"]
            + 50436 * t1708 * powerSums["20"]
            - 121296 * t855 * powerSums["23"]
            - 81804 * t855 * powerSums["41"]
            - 49208 * t858 * powerSums["32"]
            - 6 * t3086 * powerSums["23"]
            + 120 * t772 * powerSums["23"]
            + 40 * t775 * powerSums["32"]
            + 54 * t1675 * powerSums["20"]
            - 1152 * t788 * powerSums["23"]
            - 6 * t3187 * powerSums["14"]
            - 6 * t3187 * powerSums["32"]
            - 24 * t532 * powerSums["15"]
            + 120 * t898 * powerSums["14"]
            + 120 * t898 * powerSums["32"]
            + 22 * t766 * powerSums["31"]
            + 372 * t304 * powerSums["15"]
            - 1080 * t911 * powerSums["14"]
            - 1152 * t911 * powerSums["32"]
            - 2 * t2998 * powerSums["23"]
            - 6 * t3086 * powerSums["41"]
            + 40 * t530 * powerSums["23"]
            + 120 * t772 * powerSums["41"]
            + 72 * t640 * powerSums["22"]
            + 12 * t640 * powerSums["40"]
            - 1932 * t653 * powerSums["21"]
            + 526 * t2912 * powerSums["20"]
            - 2090 * t472 * powerSums["40"]
            + 6588 * t661 * powerSums["21"]
            - 1351 * t2927 * powerSums["20"]
            - 2 * t2698 * powerSums["40"]
            - 6142 * t390 * powerSums["41"]
            + 14572 * t377 * powerSums["41"]
            - 15918 * t472 * powerSums["22"]
            - 22880 * t374 * powerSums["41"]
            + 32748 * t488 * powerSums["22"]
            + 20448 * t371 * powerSums["41"]
            - 34668 * t504 * powerSums["22"]
            - 7488 * t307 * powerSums["41"]
            + 14040 * t517 * powerSums["22"]
            - t2874 * powerSums["20"]
            - 24 * t532 * powerSums["51"]
            + 32 * t415 * powerSums["40"]
            - 24 * t637 * powerSums["21"]
            + 16 * t2883 * powerSums["20"]
            - 6 * t2698 * powerSums["22"]
            - 56 * t532 * powerSums["33"]
            - 15016 * t52 * powerSums["34"]
            - 11262 * t52 * powerSums["52"]
            + 12030 * t261 * powerSums["42"]
            + 12072 * t64 * powerSums["24"]
            + 16125 * t64 * powerSums["42"]
            - 9432 * t70 * powerSums["32"]
            + 5562 * t1128 * powerSums["22"]
            - 108 * t1131 * powerSums["12"]
            + t3551
            + t3527
            + 38328 * t84 * t115
            - 18102 * t67 * t115
            + 4452 * t52 * t115
            - 480 * t33 * t115
            - 12 * t24 * t115
            + 6 * t4 * t115
            + 8964 * t84 * t5
            - 4212 * t96 * t12
            + 3474 * t70 * t21
            + 612 * t67 * t25
            - 2478 * t67 * t5
            + 3024 * t81 * t12
            + 8604 * t84 * t35
            + 51084 * t87 * t18
            - 5292 * t87 * t21
            - 1104 * t84 * t25
            + t73
            + t3477
            + t3373
            + t3451
            + t129
            + t3574
            + t3667
            + t3644
            + t3425
            + t3279
            + t3328
            + t3620
        )
        / (-7 + numberOfSamples)
        / (-6 + numberOfSamples)
        / (-5 + numberOfSamples)
        / t1413
        / t1484
        / t317
        / t3006
        / numberOfSamples
        / 4.0
    )
    return uglyD1O4E
