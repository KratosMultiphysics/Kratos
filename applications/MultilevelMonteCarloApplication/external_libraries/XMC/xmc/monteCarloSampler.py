# XMC imports
from xmc.tools import dynamicImport
from xmc.tools import splitOneListIntoTwo
from xmc.tools import mergeTwoListsIntoOne
import xmc.methodDefs_monteCarloSampler.asynchronousUpdateGlobalEstimators as mda

# External libraries
from collections import defaultdict
from itertools import chain
import warnings

from xmc.distributedEnvironmentFramework import *

class MonteCarloSampler():
    """
    This is a generic class for a multi-something Monte Carlo estimator, which
    also generates its own samples (albeit indirectly). It does not handle errors
    and tolerances (this is not an algorithm), but it does handle statistical es-
    timations, and models fitted across Monte Carlo indices.
    """

    def __init__(self,**keywordArgs):
        # Attributes
        self.indices = keywordArgs.get('indices')
        self.assemblers = keywordArgs.get('assemblers')
        self.estimatorsForAssembler = keywordArgs.get('estimatorsForAssembler')
        self.indexConstructor = dynamicImport(keywordArgs.get('indexConstructor',None))
        self.indexConstructorDictionary = keywordArgs.get('indexConstructorDictionary',None)
        self.qoiPredictor = keywordArgs.get('qoiPredictor',None)
        self.costPredictor = keywordArgs.get('costPredictor',None)
        self.estimatorsForPredictor = keywordArgs.get('estimatorsForPredictor',None)
        # TODO Why do we need costEstimatorsForPredictor?
        self.costEstimatorForPredictor = keywordArgs.get('costEstimatorsForPredictor',None)
        self.isCostUpdated = keywordArgs.get('isCostUpdated',False)
        self.errorEstimators = keywordArgs.get('errorEstimators',None)
        self.assemblersForError = keywordArgs.get('assemblersForError',None)
        # asynchronous framework settings
        self.numberBatches = keywordArgs.get('initialNumberBatches',None)
        self.batchIndices = None
        # Method Objects
        # TODO Define these below
        self.indexSet = None
        self.samples = None

    def hierarchy(self):
        extracted_hierarchy = []
        for index in self.indices:
            extracted_hierarchy.append([index.indexValue,index.sampleNumber()])
        return extracted_hierarchy

    def indexEstimation(self, coordinate, valueMethodArgs):
        return [self.indices[j].indexwiseContribution(coordinate,valueMethodArgs) for j in range(len(self.indices))]

    def indexCostEstimation(self, valueMethodArgs):
        return [self.indices[j].costMoment(valueMethodArgs) for j in range(len(self.indices))]

    def estimation(self, assemblerCoordinates=None):
        """
        For an entry of assembler, collect estimations from all indices corresponding to
        the same entry of estimatorsForAssembler. Then, pass this list to the assembleEstimation
        method of that entry of assembler.
        """

        # If nothing is specified, assemble all estimations
        if assemblerCoordinates is None:
            assemblerCoordinates = range(len(self.assemblers))

        # Extract current hierarchy from list of indices
        hierarchy = self.hierarchy()

        ## Get list of estimations for assemblers without duplicate calls

        # List of which assembler needs what estimations
        args = [self.estimatorsForAssembler[c] for c in assemblerCoordinates]
        # Convert nested list to nested tuple
        mapArg = [( (i,j), (v[0], tuple(v[1])) )
                   for i,a in enumerate(args) for j,v in enumerate(a)]
        # Create dictionary of {argument: [coord1, ...], ...}
        # such that args[coord1[0]][coord1[1]] = argument
        argMap = defaultdict(list)
        for t in mapArg:
            argMap[t[1]].append(t[0])
        # Initialise and fill list of estimations for assemblers
        estimations = [[[] for _ in a] for a in args]
        for estArgs,coords in argMap.items():
            # Compute this unique estimation
            est = self.indexEstimation(*estArgs)
            # Distribute it wherever is appeared in args
            for c in coords:
                estimations[c[0]][c[1]] = est

        # Run the corresponding estimation methods on this
        globalEstimations = []
        # Iterate over couples (coord,estimation)
        for c,e in zip(assemblerCoordinates,estimations):
            ge = self.assemblers[c].assembleEstimation(hierarchy,e)
            globalEstimations.append(ge)

        # Delete COMPSs objects
        # Flatten list of depth 2 then unpack
        delete_object(*chain.from_iterable(hierarchy),
                      *chain.from_iterable(chain.from_iterable(estimations)))

        return globalEstimations

    def errorEstimation(self, errorEstimatorCoordinates=None):
        """
        For an entry of errorEstimators, collect assembledEstimations from all assemblers
        corresponding to the same entry of assemblersForError. Then, pass this list to
        the error method of that entry of errorEstimators.
        """

        if errorEstimatorCoordinates is None:
            errorEstimatorCoordinates = range(len(self.errorEstimators))

        # Get list of estimations for assemblers without duplicate calls

        # List of which error estimator needs what estimations
        args = [self.assemblersForError[x] for x in errorEstimatorCoordinates]
        # Convert nested list to nested tuple
        mapArg = [( (i,j), v )
                   for i,a in enumerate(args) for j,v in enumerate(a)]
        # Create dictionary of {argument: [coord1, ...], ...}
        # such that args[coord1[0]][coord1[1]] = argument
        argMap = defaultdict(list)
        for t in mapArg:
            argMap[t[1]].append(t[0])
        # Initialise and fill list of global estimations for error estimation
        globalEst = [[[] for _ in a] for a in args]
        # Compute this unique estimation
        uniqueEstimations = self.estimation(list(argMap.keys()))
        # Distribute it wherever is appeard in args
        for i,coords in enumerate(argMap.values()):
            for c in coords:
                globalEst[c[0]][c[1]] = uniqueEstimations[i]

        # Compute errors
        errors = [self.errorEstimators[c].error(e)
                  for c,e in zip(errorEstimatorCoordinates,globalEst)]

        return errors

    def updateIndexSet(self, newHierarchy):
        """
        Remove entries present in self.indices that are no longer present
        in newHierarchy. Add entries present in newHierarchy that are not
        there in self.indices
        """

        newIndices,newSamples = splitOneListIntoTwo(newHierarchy)
        # Collect the positions in self.indices that are not
        # present in newHierarchy
        # TODO there must be a simpler way
        list_of_positions_to_remove = []
        for i in range(len(self.indices)):
            is_index_found = False
            for j in range(len(newIndices)):
                if(newIndices[j]==self.indices[i].indexValue):
                    is_index_found = True
                    break
            if(is_index_found is False):
                list_of_positions_to_remove.append(i)

        # Delete the positions in self.indices that are not
        # present in newHierarchy
        self._removeMonteCarloIndex(list_of_positions_to_remove)

        # Collect the entries in newHierarchy that are not
        # present in self.indices
        # TODO there must be a simpler way
        list_of_indices_to_add = []
        for i in range(len(newIndices)):
            is_index_found = False
            for j in range(len(self.indices)):
                if(newIndices[i]==self.indices[j].indexValue):
                    is_index_found = True
                    break
            if(is_index_found is False):
                list_of_indices_to_add.append(newIndices[i])

        # Pass the list of missing entries to the method that creates
        # these entries in self.indices
        self._addMonteCarloIndex(list_of_indices_to_add)

        # Rearrange self.indices to match the ordering of newHierarchy
        for i in range(len(newIndices)):
            for j in range(len(self.indices)):
                if(newIndices[i]==self.indices[j].indexValue):
                    break
            if(i==j):
                break
            else:
                self.indices[i],self.indices[j] = self.indices[j],self.indices[i]

    def _addMonteCarloIndex(self, indicesToAdd):
        """
        Creates and appends MonteCarloIndex instances corresponding to each entry of
        indicesToAdd.
        """

        assert(self.indexConstructor is not None),"Index constructor is not defined."
        for i in indicesToAdd:
            self.indexConstructorDictionary['indexValue']=i
            self.indices.append(self.indexConstructor(**self.indexConstructorDictionary))

    def _removeMonteCarloIndex(self,indicesToRemove):
        """
        Removes elements of self.indices corresponding to each entry of
        indicesToAdd.
        """

        # Warn if list is not empty
        if indicesToRemove:
            warnings.warn(('Monte Carlo Index removal requested but not expected. '
                           'I will do it anyway, but check that it is what you expect.'),
                          UserWarning)
            # Warning intended for developers, hence the category (although inaccurate)
            warnings.warn(('Future objects associated with MonteCarloIndex instances'
                           'are not removed.'),
                          PendingDeprecationWarning)
        for i in indicesToRemove:
            del self.indices[i]


    def updatePredictors(self, predictorCoordinates=None):
        """
        For an entry of qoiPredictor, retrieve all index estimations for the corresponding
        entry of estimatorsForPredictor and pass to the update method of qoiPredictor.
        If self.isCostUpdated is True, then do the same for costPredictor as well
        """

        # TODO - The implicit assumption here is that there is a one-to-one map
        # between levelwise estimations and model built upon that estimation. We
        # have not yet allowed for a mechanism similar to the assembler mechanism
        # but at the indexwise quantities. Should we construct such a mechanism?

        # If nothing is specified, update all qoiPredictors
        if predictorCoordinates is None:
            predictorCoordinates = range(len(self.qoiPredictor))

        # Extract current hierarchy from list of indices
        hierarchy = self.hierarchy()
        [indices,number_samples] = splitOneListIntoTwo(hierarchy)

        # Retrieve all index estimations corresponding to each entry of predictorCooridnates
        for i in predictorCoordinates:
            # Extract qoiEstimator coordinates and value arguements from estimatorsForPredictor
            qoi_estimator_coordinate = self.estimatorsForPredictor[i][0]
            qoi_estimator_value_arguements = self.estimatorsForPredictor[i][1]
            value_transform_function = self.estimatorsForPredictor[i][2]

            # Assemble this quantity across all indices and apply transformation
            estimations = self.indexEstimation(
                qoi_estimator_coordinate,qoi_estimator_value_arguements)
            # TODO - qoiEstimator.value returns Var(Y_l)/N_l, whereas the model is fit
            # on Var(Y_l). The following if-else is a hack to multiply by the number of
            # samples . It is here temporarily until a better mechanism can be found.
            # Refer the to-do above.
            if qoi_estimator_value_arguements[-1] is False:
                estimationsForPredictors = [value_transform_function(i) for i in estimations]
            else:
                estimationsForPredictors = [value_transform_function(index,number_samples[i])
                                                  for i,index in enumerate(estimations)]

            # Extract only the indices for which all of the multi-index components are non-zero
            data = []
            for j in range(len(indices)):
                index = indices[j]
                if all([index[i]>0 for i in range(len(index))]):
                    data.append([index,estimationsForPredictors[j]])
                else:
                    pass

            # Run update method of qoiPredictor on these index estimations
            self.qoiPredictor[i].update(data)

        # Retrieve all index cost estimations and pass them to costPredictor
        if self.isCostUpdated:
            estimationsForPredictors = self.indexCostEstimation(self.costEstimatorForPredictor)

            # Extract only the indices for which all of the multi-index components are non-zero
            data = []
            for j in range(len(indices)):
                index = indices[j]
                if all([index[i]>0 for i in range(len(index))]):
                    data.append([index,estimationsForPredictors[j]])
                else:
                    pass

            data = mergeTwoListsIntoOne(indices,estimationsForPredictors)
            self.costPredictor.update(data)

    def update(self,newHierarchy):
        """
        Add and process new samples. It first runs updateIndexSet, then calls on each indexs update, then calls updatePredictor.
        """

        # Add new indices and trim ones no longer required
        self.updateIndexSet(newHierarchy)

        for i in range(len(self.indices)):
            self.indices[i].update(newHierarchy[i])

        # synchronize estimator needed for checking convergence and updating hierarchy
        for i,index in enumerate(self.indices):
            self.indices[i].qoiEstimator = get_value_from_remote(index.qoiEstimator)
            self.indices[i].costEstimator = get_value_from_remote(index.costEstimator)

        # TODO Find a way to update predictors before sync?
        # Update models with new observations
        self.updatePredictors()


    ####################################################################################################
    ###################################### ASYNCHRONOUS FRAMEWORK ######################################
    ####################################################################################################


    def asynchronousPrepareBatches(self,newHierarchy):
        """
        Method setting-up batches. If needed, the serialize method is called, to serialize Kratos Model and Kratos Parameters. Otherwise, serialized objects are passed to new batches, to avoid serializing multiple times.
        For each batch, an index constructor dictionary is built. This way, local estimators may be computed in parallel and then update global estimators.
        """

        newIndices,newSamples = splitOneListIntoTwo(newHierarchy)
        if self.batchIndices is None: # iterationCounter = 0
            # serialze Kratos object sinto monteCarloIndex indeces instances
            # requirements: KratosSolverWrapper as SolverWrapper and concurrent_adaptive_refinement as refinement_strategy
            self.indices[0].sampler.solvers[0].serialize()
            # prepare batch indices and booleans
            self.batchIndices = [[] for _ in range (self.numberBatches)]
            self.batchesLaunched = [False for _ in range (self.numberBatches)]
            self.batchesExecutionFinished = [False for _ in range (self.numberBatches)]
            self.batchesAnalysisFinished = [False for _ in range (self.numberBatches)]
            self.batchesConvergenceFinished = [False for _ in range (self.numberBatches)]
        else: # iterationCounter > 0
            self.asynchronousUpdateBatches()
        # create monteCarloIndex instances
        for batch in range (self.numberBatches):
            if self.batchesLaunched[batch] is False:
                index = 0
                for i in newIndices:
                    self.indexConstructorDictionary['indexValue']=i
                    self.batchIndices[batch].append(self.indexConstructor(**self.indexConstructorDictionary))
                    # save needed serializations stored in indices
                    # build finer level
                    self.asynchronousSerializeBatchIndices(batch,index,solver=0)
                    # build coarser level
                    if (index > 0):
                        self.asynchronousSerializeBatchIndices(batch,index,solver=1)
                    index = index + 1


    def asynchronousSerializeBatchIndices(self,batch,index,solver):
        """
        Method passing serialized Kratos Model and Kratos Parameters and other KratosSolverWrapper members to new batches. Passing them, we avoid serializing for each batch, and we do only once for the first batch.
        """

        # model
        self.batchIndices[batch][index].sampler.solvers[solver].pickled_model = self.indices[0].sampler.solvers[0].pickled_model
        self.batchIndices[batch][index].sampler.solvers[solver].serialized_model = self.indices[0].sampler.solvers[0].serialized_model
        # project parameters
        self.batchIndices[batch][index].sampler.solvers[solver].pickled_project_parameters = self.indices[0].sampler.solvers[0].pickled_project_parameters
        self.batchIndices[batch][index].sampler.solvers[solver].serialized_project_parameters = self.indices[0].sampler.solvers[0].serialized_project_parameters
        # custom metric refinement parameters
        self.batchIndices[batch][index].sampler.solvers[solver].pickled_custom_metric_refinement_parameters = self.indices[0].sampler.solvers[0].pickled_custom_metric_refinement_parameters
        # custom remesh refinement parameters
        self.batchIndices[batch][index].sampler.solvers[solver].pickled_custom_remesh_refinement_parameters = self.indices[0].sampler.solvers[0].pickled_custom_remesh_refinement_parameters
        # booleans
        self.batchIndices[batch][index].sampler.solvers[solver].is_project_parameters_pickled = True
        self.batchIndices[batch][index].sampler.solvers[solver].is_model_pickled = True
        self.batchIndices[batch][index].sampler.solvers[solver].is_custom_settings_metric_refinement_pickled = True
        self.batchIndices[batch][index].sampler.solvers[solver].is_custom_settings_remesh_refinement_pickled = True


    def asynchronousUpdateBatches(self):
        """
        Method updating all relevant members when new batches are created.
        """

        # set here number of batches to append
        # append only if not exceeding the maximum number of iterations
        if (self.numberBatches > self.maximumNumberIterations):
            new_number_batches = 0
        else:
            new_number_batches = 1
        self.numberBatches = self.numberBatches + new_number_batches
        for _ in range (new_number_batches):
            self.batchIndices.append([])
            self.batchesLaunched.append(False)
            self.batchesExecutionFinished.append(False)
            self.batchesAnalysisFinished.append(False)
            self.batchesConvergenceFinished.append(False)


    def asynchronousInitialize(self,newHierarchy):
        """
        Method initializing new batches.
        """

        # Add new indices and trim ones no longer required
        self.updateIndexSet(newHierarchy)
        # Prepare batches
        self.asynchronousPrepareBatches(newHierarchy)


    def asynchronousLaunchEpoch(self,newHierarchy):
        """
        Method spawning new batches. It first spawn high index batch realizations, since they have a larger computational cost and normally require a larger time-to-solution.
        """

        # Run batch hierarchies
        for batch in range(self.numberBatches):
            if (self.batchesLaunched[batch] is not True):
                self.batchesLaunched[batch] = True
                for i in reversed(range(len(self.batchIndices[batch]))):
                    self.batchIndices[batch][i].update(newHierarchy[i])
                self.batchesExecutionFinished[batch] = True
                self.batchesAnalysisFinished[batch] = True


    def asynchronousFinalize(self,batch):
        """
        Method updating all global estimators with the contributions of the new batch.
        """

        # Postprocess on finished batches
        for level in range (len(self.batchIndices[batch])):
            # update global estimators
            mda.updateGlobalMomentEstimator_Task(self.indices[level].qoiEstimator,self.batchIndices[batch][level].qoiEstimator,self.indices[level].costEstimator,self.batchIndices[batch][level].costEstimator,batch)
            # delete COMPSs future objects no longer needed
            delete_object(self.batchIndices[batch][level].costEstimator,
                          *self.batchIndices[batch][level].qoiEstimator)
            # Update model coefficients for cost, bias, variance with new observations
            self.updatePredictors()

        for level in range (len(self.indices)):
            # synchronize estimator needed for checking convergence and updating hierarchy
            self.indices[level].qoiEstimator = get_value_from_remote(self.indices[level].qoiEstimator)
            self.indices[level].costEstimator = get_value_from_remote(self.indices[level].costEstimator)


    def asynchronousUpdate(self,newHierarchy):
        """
        Method called to initialize all new batches to launch, and to launch them.
        """

        self.asynchronousInitialize(newHierarchy)
        self.asynchronousLaunchEpoch(newHierarchy)
