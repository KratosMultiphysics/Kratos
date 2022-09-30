import xmc.tools as tools

from exaqute import *

# TODO When MomentEstimator is updated to specs, one level of nesting will have to be added above solver level: each element of sampleGroup is to be a list of sample components; as of now, it is just a single component.


@task(keep=True, estimators=INOUT, samples={Type: COLLECTION_IN, Depth: 3})
def updatePartialQoiEstimators_Task(estimators, samples):
    """
    Update a list of estimators with a set of samples.

    Input arguments:
    - estimators: list of StatisticalEstimator objects.
    - samples: tridimensional array (nested list) whose dimensions are (event,solver,component)

    It is assumed that estimators[i] expects the i-th component (i.e. samples[?][?][i]).
    Consequently, samples must have length len(estimators) along its last dimension.
    """

    # Iterate over estimators
    for iEst, _ in enumerate(estimators):
        # Get subset of samples for component iEst
        # as array expected by update methud: (event, solver).
        # Example:
        # samples[event][solver] = [q_0_solver, ..., q_N_solver]
        # then
        # sampleSubset[event] = [q_iEst_0, ..., q_iEst_numberOfSolvers]
        # Observe the type of each sample (float or list) depends on the chosen moment estimator
        sampleSubset = [[perSolver[iEst] for perSolver in perEvent] for perEvent in samples]
        # Pass to update method
        estimators[iEst].update(sampleSubset)


@task(keep=True, origin=IN, destination=INOUT)
def assignByIndex_Task(destination, origin, mapping):
    """
    Assign elements into a list, according to an index.

    Input arguments:
    - destination: list in which to assign elements
    - origin: list of elements to assign
    - mapping: list of integers such that origin[i] must go to destination[mapping[i]]
    """

    for i, j in enumerate(mapping):
        destination[j] = origin[i]


@task(keep=True, costEstimator=INOUT, times={Type: COLLECTION_IN, Depth: 2})
def updateCostEstimator_Task(costEstimator, times):
    """
    Update cost estimator with a set of time samples.

    Input arguments:
    - costEstimator: MomentEstimator object.
    - times: bidimensional array (nested list) of time samples, whose dimensions are
    (event,solver).
    """

    # For each event, sum times across solvers
    # Structure as expected by MomentEstimator.update: bidimensional list
    timeSums = [[sum(perEvent)] for perEvent in times]
    # Update cost estimator
    costEstimator.update(timeSums)
