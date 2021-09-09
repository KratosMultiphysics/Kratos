# Import external libraries
import itertools as it
from math import ceil
import numpy as np
import warnings

# XMC imports
from xmc.tools import instantiateObject, summation_Task
from xmc.methodDefs_monteCarloIndex import updateEstimators
from exaqute import get_value_from_remote, delete_object
import xmc.methodDefs_monteCarloIndex.updateEstimators as mdu


class MonteCarloIndex:
    """
    This is a generic class for an index of an MXMC method. It handles sample
    generation and index-specific estimators.
    """

    # TODO - this constructor should change based on the json file I/O mechanism
    # that will get set up at a future date
    def __init__(self, **keywordArgs):
        # Attributes
        self.indexValue = keywordArgs.get("indexValue")
        self.costEstimator = instantiateObject(
            keywordArgs.get("costEstimator"), **keywordArgs.get("costEstimatorInputDictionary")
        )
        # TODO 'qoiEstimator' should be renamed
        self.qoiEstimator = []
        # qoi estimators
        qoi_estimator_module = keywordArgs.get("qoiEstimator")
        qoi_estimator_module_args = keywordArgs.get("qoiEstimatorInputDictionary")
        for i, d in enumerate(qoi_estimator_module_args):
            self.qoiEstimator.append(instantiateObject(qoi_estimator_module[i], **d))

        sampler_input_dict = keywordArgs.get("samplerInputDictionary")
        sampler_input_dict["solverWrapperIndices"] = self.solverIndices()
        self.sampler = instantiateObject(keywordArgs.get("sampler"), **sampler_input_dict)
        # TODO remove. Samples must be stored in self.data.
        self.samples = keywordArgs.get("samples", None)
        self.areSamplesStored = keywordArgs.get("areSamplesStored", False)
        self.areSamplesRecycled = keywordArgs.get("areSamplesRecycled", True)
        self.eventGroupSize = keywordArgs.get("eventGroupSize", None)

    def sampleNumber(self):
        if self.areSamplesStored:
            return len(self.samples)
        else:
            return self.qoiEstimator[
                0
            ].sampleNumber()  # necessary for minimizing synchronization points

    # TODO - potentially better name here
    def indexwiseContribution(self, qoiEstimatorCoordinate, qoiEstimatorValueArguements):
        return self.qoiEstimator[qoiEstimatorCoordinate].value(*qoiEstimatorValueArguements)

    def costMoment(self, costEstimatorValueArguments):
        """
        Returns the requested moment estimation of the cost.
        The input arguments are passed to the cost estimator's relevant method.
        """

        return self.costEstimator.value(*costEstimatorValueArguments)

    def newSample(self):
        """
        Draw  a single sample, i.e. set of correlated outputs from the same event.

        See the documentation of SampleGenerator.generate for the data structure.
        """

        sample = self.sampler.generate()
        return sample

    def indexSetDimension(self):
        """
        Returns the dimension of the Monte Carlo index set.
        """

        return len(self.indexValue)

    def solverIndices(self):
        """
        Compute list of indices of correlated solvers.
        """

        # TODO What follows must be reviewed and probably revised.

        diff_tuple = []
        hypercube_vertices = 2 ** self.indexSetDimension()
        for _ in range(self.indexSetDimension()):
            diff_tuple.append([0, 1])
        diff_list = list(it.product(*diff_tuple))
        for i in range(len(diff_list)):
            diff_list[i] = list(diff_list[i])

        main_index = [self.indexValue] * hypercube_vertices
        return [
            [main_index[i][j] - diff_list[i][j] for j in range(self.indexSetDimension())]
            for i in range(hypercube_vertices)
        ]

    def _addSamples(self, newSamples):
        """
        Add newsamples to stored samples. This is an internal method with no check on whether sample storage is enabled; be careful when using it.
        """

        # TODO replace self.samples with self.data
        self.samples.append(newSamples)

    def update(self, newIndexAndSampleNumbers):
        """
        Update the Monte Carlo index to a new number of samples. First,
        decide the number of new samples to be generated, then generate
        each sample and the associated cost with newSample and pass them
        to the update methods of qoiEstimator and costEstimator.
        """

        # Compute number of new samples based on areSamplesRecycled
        # TODO Minimum number of new samples hard coded here to 6, since
        # program not planned for moment > 4 estimation. Accept/infer this
        # value later based on max(qoiEstimator[i].order)
        if self.areSamplesRecycled is True:
            number_new_samples = max(5, newIndexAndSampleNumbers[1] - self.sampleNumber())
        else:
            number_new_samples = newIndexAndSampleNumbers[1]
            for i in range(len(self.qoiEstimator)):
                self.qoiEstimator[i].reset()
            self.costEstimator.reset()

        ### Drawing samples

        # Generate the required number of correlated samples
        # and estimate cost
        samples = []
        times = []
        for _ in range(number_new_samples):
            # Generate a new sample (list of solver outputs for one random event)
            # See the documentation of method newSample for the data structure.
            new_sample, new_time = self.newSample()
            # append to corresponding list
            samples.append(new_sample)
            times.append(new_time)

        # Note on data structure of samples
        # ---------------------------------

        # WARNING: This is adapted from the documentation of SampleGenerator and
        # may be outdated. It is recommended that you use the documentation instead.

        # Let us note S = samples, S_j = samples[j], etc.. S is a nested lists of depth > 2.
        # len(S) == number_new_samples; len(S_j) == MonteCarloIndex.numberOfSolvers().
        # Elements of S_jk may be future objects.

        # If solver outputs are not split, (i.e. not MonteCarloIndex.areSamplesSplit()),
        # len(S_jk) == MonteCarloIndex.sampleDimension().
        # Example, for any event i. If MonteCarloIndex.numberOfSolvers() == 2 and
        # SampleGenerator.sampleDimension() == 3:
        # S_i1 == [s1_o1, s1_o2, s1_o3] (== [future, future, future] if parallel)
        # (s: solver, o: output); idem for S_i2.

        # If solver outputs are split, (i.e. MonteCarloIndex.areSamplesSplit()),
        # S_ij is of length len(MonteCarloIndex.sampleSplitSizes())
        # and S_ijk is of length MonteCarloIndex.sampleSplitSizes()[k].
        # Example, for any event i. If MonteCarloIndex.numberOfSolvers() == 2 and
        # MonteCarloIndex.sampleSplitSizes() == [3, 2]:
        # S_i1 == [ S_(1,1), S_(1,2) ] == [ [s1_o1, s1_o2, s1_o3], [s1_o4, s1_o5] ]
        #      == [future, future] if parallel
        # (s: solver, o: output); idem for solver S_i2.

        ### Estimator update

        # Get the multi-indices to the subsets of estimators, samples and times
        # that will be used to update the estimators
        indexEst, indexSamp, indexTime = self._updateSubsets(number_new_samples)
        # Iterator over the 'solver' dimension of the sample and time arrays
        # TODO Idea: include this in multi-indices returned by _updateSubsets?
        solverRg = range(self.numberOfSolvers())

        ## Case: solver outputs are not split

        if not self.areSamplesSplit():
            # Iterate over subsets of estimators
            # For the moment, we actually iterate over estimators (subsets are singletons)
            for g, iE in enumerate(indexEst):
                # Update for self.qoiEstimators
                # Iterate over subsets of events
                for indexS in indexSamp[g]:
                    # Assemble subset of samples
                    # TODO When MomentEstimator is updated to specs, one level of nesting will
                    # have to be added above solver level: each element of sampleGroup is to
                    # be a list of sample components; as of now, it is just a single component.
                    sampleGroup = [[samples[i][j][k] for j in solverRg] for i, k in indexS]
                    # Update estimator with sample subset
                    self.qoiEstimator[iE].update(sampleGroup)

            # Update self.costEstimator
            # Iterate over subsets of events
            for indexT in indexTime:
                # Assemble subset of samples
                # in bidimensional array expected by MomentEstimator.update
                # TODO There must be a better way to handle this
                timeGroup = [[summation_Task(*times[i])] for i in indexT]
                # Update estimator with sample subset
                self.costEstimator.update(timeGroup)

            return  # What follows is executed only if self.areSamplesSplit()

        ## Case: solver outputs are split
        # Iterate over splits of solver outputs
        for g, iE in enumerate(indexEst):
            # Assemble subset of estimators for current split
            estimatorGroup = [self.qoiEstimator[ie] for ie in iE]
            # Update for self.qoiEstimators
            # Iterate over subsets of events
            for indexS in indexSamp[g]:
                # Assemble subset of samples
                sampleGroup = [[samples[i][j][k] for j in solverRg] for i, k in indexS]
                # Update subset of estimators
                mdu.updatePartialQoiEstimators_Task(estimatorGroup, sampleGroup)
                # Delete future objects no longer needed
                # Pass an iterable over *future* objects (e.g. no list of futures)
                # therefore, flatten list of depth 2
                delete_object(*it.chain.from_iterable(sampleGroup))
            # Re-assign estimators from updated subset
            mdu.assignByIndex_Task(self.qoiEstimator, estimatorGroup, iE)
            # Delete future objects no longer needed
            delete_object(*estimatorGroup)

        # Update self.costEstimator
        # Iterate over subsets of events
        for indexT in indexTime:
            # Assemble subset of samples
            # as list expect by updateCostEstimator_Task
            timeGroup = [times[i] for i in indexT]
            # Update cost estimator with sample subset
            # TODO Unnecessary, we should use MomentEstimator.update
            mdu.updateCostEstimator_Task(self.costEstimator, timeGroup)
            # delete future objects no longer needed
            delete_object(*it.chain.from_iterable(timeGroup))

    def _updateSubsets(self, numberOfSamples: int):
        """
        Returns the multi-indices of the subsets – of estimators, samples and times – used for update.

        Input arguments:
        - number of samples: int, number of random events = length of list of samples = number of calls to the sampler.

        Output arguments:
        - multi-index of estimator subsets: nested list
        - multi-index of sample subsets: nested list
        - multi-index of time (cost) subsets: nested list


        We update each estimator by passing it samples. These samples are passed by groups of fixed size (self.eventGroupSize), corresponding to a subset of events.
        Each estimator needs only some components. Moreover, the outputs from the solvers may be split in sublists.
        Therefore, the sets of samples and time received from the sample generator, as well as the list of estimator, need to be re-arranged into subsets suitable for this update process.
        This method generates the multi-indices used to extract the necessary subsets from these multi-dimensional arrays.
        We use below the notation A_ij := A_(i,j) := A[i][j] and (A_ij) the multi-dimensional array of these components for all acceptable values of i and j.

        ## Subsets of times (costs)

        The array of times T := (T_ij) is bidimensional and its dimensions are (event, solver): T_ij is the computational time returned by the j-th solver for the i-th random event.
        We want to create from T a tridimensional array gT = (gT_ijk) of dimensions (subset,event,solver) such that gT_i is the i-th subset passed to the cost estimator.
        This method creates (among others) the multi-index iT such that gT = T_iT; iT is bidimensional.
        The notation gT = T_iT is to be understood as: iT_ijk is an index of appropriate dimension for T (i.e. 2) and gT_ijk = T_(iT_ijk) for all acceptable values of i, j and k.
        E.g. if gT_(4,1,3) = T_(2,18), then iT_(4,1,3) = (2,18).

        ## Subsets of samples and estimators

        # Case: solver outputs are not split

        If solver outputs are not split, the array of samples S := (S_ijk) is tridimensional.
        Its dimensions are (event,solver,component): S_ijk is the k-th component (or output, or random variable) of the j-th solver for the i-th random event.
        We wish to create from the array of estimators E a unidimensional array (gE_i), and from S a tetradimensional array gS = (gS_ijklm),
        so that gS_ij is the j-th subset to be passed to gE_i. The dimensions of gS are (subset,component,event,solver), and those of gE are (component).
        In this case, gE = E; not so if solver outputs are split (see below).
        This method creates the multi-indices iE, iS such that gS = S_iS and gE = E_iE – in addition to the multi-index of the time subset (see above).
        Since the last dimension (i.e. solver) is unchanged from S to gS, iS is only tridimensional. Obviously, iE is unidimensional.

        # Case: solver outputs are split

        If solver outputs are split, S := (S_ijkl) is quadridimensional, and its dimensions are (event,solver,split,component).
        However, in this case the sample subsets are not passed to the update method of a StatisticalEstimator, but to another method (updatePartialQoiEstimators_Task).
        Unlike StatisticalEstimator.update, who expects an array of dimensions (event,solver), this method expects an array of dimensions (event,split,solver,component).
        Therefore, we can ignore the component dimension in arranging the subsets of samples: we consider S := (S_ijk) of dimensions (event, solver, split).
        Unlike the previous case, gE is now a bidimensional array of dimensions (splits,component): gE_ij is the estimator receiving the j-th component of the i-th split.
        We define iE and iS as previously: gE = E_iE and gS = S_iS.
        """

        # TODO Idea: make multi-indices nested tuples instead of nested lists?

        eventGpSz = self.eventGroupSize
        if eventGpSz is None:
            # This is the default: update with all samples at once
            eventGpSz = numberOfSamples

        # The part that depends on whether solver outputs are split
        if self.areSamplesSplit():
            splitSz = self.sampleSplitSizes()
            # Estimator subset
            # Split the list of estimator the same way as solver outputs
            s = np.cumsum([0, *splitSz])  # indices at which solver outputs are split
            estimatorSubset = [list(range(s[i], s[i + 1])) for i, _ in enumerate(s[:-1])]
            # Iterator for third dimension of sample array: split, here.
            splitRg = range(len(splitSz))
        else:
            # Estimator subset is the whole set, so multi-index is the identity
            estimatorSubset = list(range(len(self.qoiEstimator)))
            # Iterator for third dimension of sample array: component, here.
            splitRg = range(self.sampleDimension())

        # Sample and time subsets
        numberOfEventGroups = ceil(numberOfSamples / eventGpSz)
        # List event indices delimiting the groups for update
        # The min makes sure the number of event is not exceeded
        # (i.e. last group may be smaller than the others)
        eventGroups = [
            min(i * eventGpSz, numberOfSamples) for i in range(numberOfEventGroups + 1)
        ]
        # Initialise empty nested lists for multi-indices
        sampleSubset = [[[] for _ in range(numberOfEventGroups)] for _ in splitRg]
        timeSubset = [[] for _ in range(numberOfEventGroups)]
        for b in range(numberOfEventGroups):
            # Range of event indices for current subset
            eventRg = range(eventGroups[b], min(eventGroups[b + 1], numberOfSamples))
            # Iterate over estimator subsets
            for k in splitRg:
                # For estimator subset (k,b),
                # store event index and third-dimension (split or component) index
                sampleSubset[k][b] = [(i, k) for i in eventRg]
            # For time, we know the estimator and the component,
            # just store event indices
            timeSubset[b] = list(eventRg)

        return estimatorSubset, sampleSubset, timeSubset

    def numberOfSolvers(self):
        """
        Returns the numberofsolvers in the sample generator<
        """

        return len(self.sampler.solvers)

    def sampleDimension(self, *args) -> int:
        """
        Returns the number of outputs from the sampler.
        This is the number of outputs from a single solver of the sampler,
        with no information regarding splitting, if any.

        Input arguments: same as SampleGenerator.sampleDimension

        Output argument:
        - sample dimension: integer. See SampleGenerator.sampleDimension
        """

        return self.sampler.sampleDimension(*args)

    def sampleSplitSizes(self, *args):
        """
        Returns the lengths of (sub-)lists into which samples are split, if they are; returns None if they are not.

        Input arguments: same as SampleGenerator.samplesSplitSizes.

        Output argument: list of integers, or None. If list, it is [len(subSample) for subSample in sample], where sample is the output from a single solver (after processing).
        """

        return self.sampler.sampleSplitSizes(*args)

    def areSamplesSplit(self, *args) -> bool:
        """
        Returns True if samples are split into future (sub-)lists, False otherwise.
        Accepts the same input arguments as SampleGenerator.areSamplesSplit.

        Input arguments: same as SampleGenerator.areSamplesSplit.

        Output argument: boolean, True is samples are split and False otherwise.
        """

        return self.sampler.areSamplesSplit(*args)
