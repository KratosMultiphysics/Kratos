import warnings
import numpy as np
from xmc.distributedEnvironmentFramework import (
    ExaquteTask,
    Type,
    Depth,
    COLLECTION_IN,
    INOUT,
)
from xmc.tools import flatten
from .types import (
    PowerSumsDict,
    SampleArray,
    PowerSumsDictUL,
    CombinedSampleArray,
)
from typing import Union


def updatedPowerSums(powerSums: PowerSumsDict, samples: SampleArray) -> PowerSumsDict:
    """
    Increments power sums from an array of samples.

    Supports multi-valued random variables, and Monte Carlo index sets of dimension up to 1.

    Input arguments:
    - powerSums: dictionary following the format of MomentEstimator._powerSums
    - samples: array of samples (see below).

    Output arguments:
    - powerSums: same as input, with updated values.

    The array of samples is expected to be an array of dimension at least three, viz.
    1. random event;
    2. Monte Carlo index (e.g. solver, model or level);
    3. (and beyond) random variable component.

    Any dimension beyond the third will be collapsed into this latter one.
    """
    # Ensure proper format
    # Convert to NumPy multi-dimensional arrays
    samples = np.array(samples)
    # Ensure that this is a tridimensional array
    if len(samples.shape) > 3:
        # Assumption: components of the multi-valued random variable were split into sub-lists.
        # Unpacked these sublist by reshaping:
        # all dimensions beyond the 2nd one are collapsed into one.
        warnings.warn(
            (
                "The samples are nested too deep. I assume that you misused the split-update"
                " mechanism, and will try to fix the data structure. However, this behaviour"
                " is only supported until 2020-11. Consider using the normal update mechanism"
                " or fixing your solver wrapper."
            ),
            FutureWarning,
        )
        samples = samples.reshape((*samples.shape[:2], -1))
    elif len(samples.shape) < 3:
        # This should not happen
        raise ValueError(
            f"Input argument has {len(samples.shape)} dimensions, "
            "whereas I expected at least 3."
        )
    # Proper format has been ensured.
    # Now check the dimension of the Monte Carlo index set
    if samples.shape[1] == 2:
        # Operator of differences for an index set of dimension 1
        diffOperator = np.array([[1, 1], [1, -1]])
        # Apply the operator to the second dimension of the array of samples
        # then restore the order of dimensions
        samples = np.tensordot(diffOperator, samples, axes=([1], [1])).transpose(1, 0, 2)
    # Iterate over keys to compute each corresponding power sum
    for key in powerSums.keys():
        # Deduce array of powers from key (e.g. '21' -> [2,1]) then broadcast it
        # along the second dimension of a NumPY array of dimension 3
        powers = np.array(list(map(int, key))).reshape((1, -1, 1))
        # Increment the power sum.
        # Elementwise exponentiation and multiplication across dimension 2 (levels)
        # then summation across dimension 1 (random events).
        powerSums[key] += np.sum(np.prod(samples ** powers, axis=1), axis=0)
    return powerSums


@ExaquteTask(
    samples={Type: COLLECTION_IN, Depth: 3},
    powerSums=INOUT,
)
def updatedPowerSums_Task(powerSums: PowerSumsDict, samples: SampleArray) -> PowerSumsDict:
    return updatedPowerSums(powerSums, samples)


def addPowerSumsAndCounter(
    psDict: Union[PowerSumsDict, PowerSumsDictUL],
    counter: int,
    psArray: CombinedSampleArray,
    keyOrder: tuple = None,
) -> (Union[PowerSumsDict, PowerSumsDictUL], int):
    """
    Increment existing dictionary of power sums and sample counter for an array of such. This
    function is meant as a defintion for MultiCombinedPowerSums._powerSumsUpdater.

    Input argument:
    - psDict: dictionary of power sums to be incremented; in the format of
    MultiCombinedPowerSums._powerSums
    - psArray: multidimensional array of power sums and sample count to add; in the format
    expected by MultiCombinedMomentEstimator.update.
    - counter: sample counter to be incremented.
    - keyOrder: (optional) keys of entries of power sums in psDict, in their order in psArray.
    Default value: ('1', '2', ..., 'p'), with 'p' the number of distinct power sums in psArray.

    Output arguments:
    - psDict: same format as input (depends on dimension of the MC index set)
    - counter: same format as input

    See the documentation of MultiCombinedPowerSums for more details.
    """
    # First consider the consider the case of samples from several MC levels
    if "upper" in psDict.keys():
        # If the index set dimension is 1, then psDict is of the form
        # {'upper': psDictUpper, 'lower': psDictLower}
        # So we recurse over each sub-dictionary
        #
        # First, check that we have the expected keys
        if sorted(psDict.keys()) != sorted(("lower", "upper")):
            raise ValueError(
                "Expected the dictionary of power sums to have keys ('upper', 'lower'). "
                f"Found {tuple(psDict.keys())} instead."
            )
        # Get samples for upper level only (first element along second axis)
        psArrayUpper = [[oneEvent[0]] for oneEvent in psArray]
        # Update power sums for upper level by recursion
        psDict["upper"], counterUpper = addPowerSumsAndCounter(
            psDict["upper"], counter, psArrayUpper, keyOrder
        )
        # Now do the same for lower level
        psArrayLower = [[oneEvent[1]] for oneEvent in psArray]
        psDict["lower"], counterLower = addPowerSumsAndCounter(
            psDict["lower"], counter, psArrayLower, keyOrder
        )
        # Check that the new counters are equal and return updated arguments
        # Special case: if the lower level is fictitious (dummy sample), then
        # its counter is expected to be None and a warning has been issue (see below).
        # This case will removed onced this workaround is not supported any more.
        is_lower_dummy = counterLower is None
        if not (counterUpper == counterLower or is_lower_dummy):
            raise ValueError(
                "Expected upper and lower levels to have equal sample count. "
                f"Received {counterUpper-counter} and {counterLower-counter}, respectively."
            )
        return psDict, counterUpper
    #
    # Code below is for index set of dimension 0
    # This is reached either because this is a single-level moment estimation,
    # or because this is a recursion entered from the conditional statement above
    #
    # Try to detect if this is actually a dummy sample
    is_at_least_4d = hasattr(psArray[0][0][0], "__iter__")
    has_nonzero_value = any(flatten(psArray))
    if not is_at_least_4d and not has_nonzero_value:
        # This is a sample full of zeros, with the wrong shape.
        # We are into a recursion. Warn about future deprecation
        warnings.warn(
            (
                "Dummy samples simulating a solver level are not needed here. "
                "They will be deprecated in the near future."
            ),
            FutureWarning,
        )
        # Do nothing and return power sums as they are.
        # Return counter as None, to let caller know that this is a dummy sample.
        # Note: the caller is this function, from the conditional statement above.
        return psDict, None
    #
    # Assign default value of keyOrder, if undefined
    if keyOrder is None:
        # Default key order
        # Let us not assume the maximal order of power sums
        # We substract 1 because of the counter
        keyOrder = tuple(str(i + 1) for i in range(len(psArray[0][0][0]) - 1))
    # Check that all keys exist
    if sorted(keyOrder) != sorted(psDict.keys()):
        raise ValueError(
            "Failed to match keys of new power sums and existing ones: "
            f"{sorted(keyOrder)} versus {sorted(psDict.keys())}."
        )
    # Let's reformat psArray
    # Collapse 'solver' (i.e. second) axis, since its length is 1
    psArray = [oneEvent[0] for oneEvent in psArray]
    # Increment counter and remove count from psArray
    for oneEvent in psArray:
        # We remove the counter from every component
        for oneComponent in oneEvent:
            increment = oneComponent.pop(-1)
        # However, we increment the counter only once
        counter += increment
    # Now we expect to have a homogenous tri-dimensional array: event, component, power
    # There may be an extra axis, e.g. because of an unwanted split of solver output
    # We will convert it to a NumPy array and collapse any axis beyond the third
    psArray = np.array(psArray)
    if len(psArray.shape) > 3:
        psArray = psArray.reshape((*psArray.shape[:2], -1))
    # Now we surely have a tri-dimensional array of axes: event, component, power
    # Sum over 'events' (i.e. first) axis
    psArray = np.sum(psArray, axis=0)
    # We permute the axes to be (power, component) instead of (component, power)
    psArray = psArray.transpose()
    # Increment power sums for each power
    for i, key in enumerate(keyOrder):
        # This is actually a vector operation, using NumPy's vector addition
        psDict[key] += psArray[i]
    return psDict, counter


# The depth value is necessary here
@ExaquteTask(psDict=INOUT, psArray={Type: COLLECTION_IN, Depth: 4})
def addPowerSumsAndCounter_Task(
    psDict: Union[PowerSumsDict, PowerSumsDictUL],
    counter: int,
    psArray: CombinedSampleArray,
    keyOrder: tuple = None,
) -> (Union[PowerSumsDict, PowerSumsDictUL], int):
    return addPowerSumsAndCounter(psDict, counter, psArray, keyOrder)
