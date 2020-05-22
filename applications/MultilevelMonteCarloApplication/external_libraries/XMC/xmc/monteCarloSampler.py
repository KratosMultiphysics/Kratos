# XMC imports
from xmc.tools import dynamicImport
from xmc.tools import getUnionAndMap
from xmc.tools import splitOneListIntoTwo
from xmc.tools import mergeTwoListsIntoOne
from xmc.tools import unpackedList
import xmc.methodDefs_monteCarloSampler.asynchronousUpdate as mcsua

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
        extracted_hierarchy = self.hierarchy()

        # Construct union of estimators required, as well as corresponding map
        estimator_union, estimator_union_map = getUnionAndMap([self.estimatorsForAssembler[x] for x in assemblerCoordinates])

        # Retrieve all index estimations corresponding to each entry of estimator_union
        index_estimation_union = []
        for coordinate_and_args in estimator_union:
            qoi_estimator_coordinate = coordinate_and_args[0]
            qoi_estimator_value_arguements = coordinate_and_args[1]

            # Assemble this quantity across indices
            index_estimation_for_assembler = self.indexEstimation(qoi_estimator_coordinate,qoi_estimator_value_arguements)
            index_estimation_union.append(index_estimation_for_assembler)

        # Remap this union back onto original estimatorsForAssembler list
        index_estimation_lists = []
        for map_list in estimator_union_map:
            index_estimation_lists.append([index_estimation_union[x] for x in map_list])

        # Run the corresponding estimation methods on this
        estimations = []
        for i in range(len(assemblerCoordinates)):
            estimations.append(self.assemblers[assemblerCoordinates[i]].assembleEstimation( extracted_hierarchy,*unpackedList(index_estimation_lists[i])))
        return estimations

    def errorEstimation(self, errorEstimatorCoordinates=None):
        """
        For an entry of errorEstimators, collect assembledEstimations from all assemblers
        corresponding to the same entry of assemblersForError. Then, pass this list to
        the error method of that entry of errorEstimators.
        """

        if errorEstimatorCoordinates is None:
            errorEstimatorCoordinates = range(len(self.errorEstimators))

        # Construct union of assemblers required, as well as corresponding assembler map
        assembler_union, assembler_union_map = getUnionAndMap([self.assemblersForError[x] for x in errorEstimatorCoordinates])

        # Obtain estimations corresponding to assembler_union
        assembled_estimations_before_map = self.estimation(assembler_union)

        # Remap assembled_estimations_before_map onto original assembler data structure
        assembled_estimation_lists = []
        for map_list in assembler_union_map:
            assembled_estimation_lists.append([assembled_estimations_before_map[x] for x in map_list])

        # Compute errors
        errors = []
        for i in range(len(errorEstimatorCoordinates)):
            coord = errorEstimatorCoordinates[i]
            assembled_estimations = assembled_estimation_lists[i]
            errors.append(self.errorEstimators[coord].error(*unpackedList(assembled_estimations)))
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
        for i in range(len(list_of_positions_to_remove)):
            del self.indices[list_of_positions_to_remove[i]]

        # Collect the entries in newHierarchy that are not
        # present in self.indices
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
        self.newMonteCarloIndex(list_of_indices_to_add)

        # Rearrange self.indices to match the ordering of newHierarchy
        for i in range(len(newIndices)):
            for j in range(len(self.indices)):
                if(newIndices[i]==self.indices[j].indexValue):
                    break
            if(i==j):
                break
            else:
                self.indices[i],self.indices[j] = self.indices[j],self.indices[i]

    def newMonteCarloIndex(self, listOfIndicesToAdd):
        """
        Create and return MonteCarloIndex instances corresponding to each entry of
        listOfIndicesToAdd.
        """
        assert(self.indexConstructor is not None),"Index constructor is not defined."
        new_indices = []
        for i in listOfIndicesToAdd:
            self.indexConstructorDictionary['indexValue']=i
            self.indices.append(self.indexConstructor(**self.indexConstructorDictionary))

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
        extracted_hierarchy = self.hierarchy()
        [indices,number_samples] = splitOneListIntoTwo(extracted_hierarchy)

        # Retrieve all index estimations corresponding to each entry of predictorCooridnates
        for i in predictorCoordinates:
            # Extract qoiEstimator coordinates and value arguements from estimatorsForPredictor
            qoi_estimator_coordinate = self.estimatorsForPredictor[i][0]
            qoi_estimator_value_arguements = self.estimatorsForPredictor[i][1]
            value_transform_function = self.estimatorsForPredictor[i][2]

            # Assemble this quantity across all indices and apply transformation
            index_estimations_before_transform = self.indexEstimation(qoi_estimator_coordinate,qoi_estimator_value_arguements)
            # TODO - qoiEstimator.value returns Var(Y_l)/N_l, whereas the model is fit
            # on Var(Y_l). The following if-else is a hack to multiply by the number of
            # samples . It is here temporarily until a better mechanism can be found.
            # Refer the to-do above.
            if qoi_estimator_value_arguements[-1] is False:
                index_estimation_for_predictor = [value_transform_function(index_estimations_before_transform[i])
                                for i in range(len(index_estimations_before_transform))]
            else:
                index_estimation_for_predictor = [value_transform_function(index_estimations_before_transform[i],number_samples[i])
                                for i in range(len(index_estimations_before_transform))]

            # Extract only the indices for which all of the multi-index components are non-zero
            data = []
            for j in range(len(indices)):
                index = indices[j]
                if all([index[i]>0 for i in range(len(index))]):
                    data.append([index,index_estimation_for_predictor[j]])
                else:
                    pass

            # Run update method of qoiPredictor on these index estimations
            self.qoiPredictor[i].update(*unpackedList(data))

        # Retrieve all index cost estimations and pass them to costPredictor
        if self.isCostUpdated:
            index_estimation_for_predictor = self.indexCostEstimation(self.costEstimatorForPredictor)

            # Extract only the indices for which all of the multi-index components are non-zero
            data = []
            for j in range(len(indices)):
                index = indices[j]
                if all([index[i]>0 for i in range(len(index))]):
                    data.append([index,index_estimation_for_predictor[j]])
                else:
                    pass

            data = mergeTwoListsIntoOne(indices,index_estimation_for_predictor)
            self.costPredictor.update(*unpackedList(data))

    def update(self,newHierarchy):
        """
        Add and process new samples. It first runs updateIndexSet, then calls on each indexs update, then calls updatePredictor.
        """

        # Add new indices and trim ones no longer required
        self.updateIndexSet(newHierarchy)

        for i in range(len(self.indices)):
                self.indices[i].update(newHierarchy[i])

        # Update model coefficients for cost, bias, variance with new
        # observations
        self.updatePredictors()


    ####################################################################################################
    ###################################### ASYNCHRONOUS FRAMEWORK ######################################
    ####################################################################################################


    def asynchronousPrepareBatches(self,newHierarchy):
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
        # set here number of batches to append
        # append only if not exceeding the maximum number of iterations
        if (self.numberBatches > self.maximumNumberIterations):
            new_number_batches = 0
        else:
            new_number_batches = 1
        self.numberBatches = self.numberBatches + new_number_batches
        for new_batch in range (new_number_batches):
            self.batchIndices.append([])
            self.batchesLaunched.append(False)
            self.batchesExecutionFinished.append(False)
            self.batchesAnalysisFinished.append(False)
            self.batchesConvergenceFinished.append(False)

    def asynchronousInitialize(self,newHierarchy):
        # Add new indices and trim ones no longer required
        self.updateIndexSet(newHierarchy)
        # Prepare batches
        self.asynchronousPrepareBatches(newHierarchy)

    def asynchronousLaunchEpoch(self,newHierarchy):
        # Run batch hierarchies
        for batch in range(self.numberBatches):
            if (self.batchesLaunched[batch] is not True):
                self.batchesLaunched[batch] = True
                for i in reversed(range(len(self.batchIndices[batch]))):
                    self.batchIndices[batch][i].update(newHierarchy[i])
                self.batchesExecutionFinished[batch] = True
                self.batchesAnalysisFinished[batch] = True

    def asynchronousFinalize(self,batch):
        # Postprocess on finished batches
        for level in range (len(self.batchIndices[batch])):
            # Update global qoi estimators
            for qoi_index in range (self.batchIndices[batch][level].numberOfSamplerOutputs()):
                if (batch == 0):
                    self.indices[level].qoiEstimator[qoi_index].powerSums = self.batchIndices[batch][level].qoiEstimator[qoi_index].powerSums
                    self.indices[level].qoiEstimator[qoi_index]._sampleCounter = self.batchIndices[batch][level].qoiEstimator[qoi_index]._sampleCounter
                else:
                    # asynchronous MC
                    if (self.indices[level].qoiEstimator[qoi_index].order == 1 and self.indices[level].qoiEstimator[qoi_index].indexSetDimension == 0):
                        globalIndexPowerSum1,globalIndexPowerSum2 = mcsua.updateGlobalMonteCarloIndexOrder2Dimension0_Task(\
                            self.indices[level].qoiEstimator[qoi_index].powerSums[0][0],\
                            self.indices[level].qoiEstimator[qoi_index].powerSums[1][0],\
                            self.batchIndices[batch][level].qoiEstimator[qoi_index].powerSums[0][0],\
                            self.batchIndices[batch][level].qoiEstimator[qoi_index].powerSums[1][0])
                        self.indices[level].qoiEstimator[qoi_index].powerSums[0][0] = globalIndexPowerSum1
                        self.indices[level].qoiEstimator[qoi_index].powerSums[1][0] = globalIndexPowerSum2
                        del(globalIndexPowerSum1,globalIndexPowerSum2)
                    # asynchronous MLMC
                    elif (self.indices[level].qoiEstimator[qoi_index].order == 1 and self.indices[level].qoiEstimator[qoi_index].indexSetDimension == 1):
                        globalIndexPowerSum10,globalIndexPowerSum01,globalIndexPowerSum20,globalIndexPowerSum11,globalIndexPowerSum02 = mcsua.updateGlobalMonteCarloIndexOrder2Dimension1_Task(\
                            self.indices[level].qoiEstimator[qoi_index].powerSums[0][0],\
                            self.indices[level].qoiEstimator[qoi_index].powerSums[0][1],\
                            self.indices[level].qoiEstimator[qoi_index].powerSums[1][0],\
                            self.indices[level].qoiEstimator[qoi_index].powerSums[1][1],\
                            self.indices[level].qoiEstimator[qoi_index].powerSums[1][2],\
                            self.batchIndices[batch][level].qoiEstimator[qoi_index].powerSums[0][0],\
                            self.batchIndices[batch][level].qoiEstimator[qoi_index].powerSums[0][1],\
                            self.batchIndices[batch][level].qoiEstimator[qoi_index].powerSums[1][0],\
                            self.batchIndices[batch][level].qoiEstimator[qoi_index].powerSums[1][1],\
                            self.batchIndices[batch][level].qoiEstimator[qoi_index].powerSums[1][2])
                        # update
                        self.indices[level].qoiEstimator[qoi_index].powerSums[0][0] = globalIndexPowerSum10
                        self.indices[level].qoiEstimator[qoi_index].powerSums[0][1] = globalIndexPowerSum01
                        self.indices[level].qoiEstimator[qoi_index].powerSums[1][0] = globalIndexPowerSum20
                        self.indices[level].qoiEstimator[qoi_index].powerSums[1][1] = globalIndexPowerSum11
                        self.indices[level].qoiEstimator[qoi_index].powerSums[1][2] = globalIndexPowerSum02
                        del(globalIndexPowerSum10,globalIndexPowerSum01,globalIndexPowerSum20,globalIndexPowerSum11,globalIndexPowerSum02)
                    else:
                        raise Exception ("Order or dimension not supported. Add the required task into xmc/methodDefs_monteCarloSampler/asynchronousUpdate.py")
                    self.indices[level].qoiEstimator[qoi_index]._sampleCounter = self.indices[level].qoiEstimator[qoi_index]._sampleCounter + self.batchIndices[batch][level].qoiEstimator[qoi_index]._sampleCounter

            # Update global cost estimator
            if (batch == 0):
                self.indices[level].costEstimator.powerSums = self.batchIndices[batch][level].costEstimator.powerSums
                self.indices[level].costEstimator._sampleCounter = self.batchIndices[batch][level].costEstimator._sampleCounter
            else:
                for order in range (len(self.indices[level].costEstimator.powerSums)):
                        if (len(self.indices[level].costEstimator.powerSums[order]) == 2): # order power sum is two
                            globalIndexPowerSum1,globalIndexPowerSum2 = mcsua.updateGlobalMonteCarloIndexOrder2Dimension0_Task(self.indices[level].costEstimator.powerSums[order][0],self.indices[level].costEstimator.powerSums[order][1],self.batchIndices[batch][level].costEstimator.powerSums[order][0],self.batchIndices[batch][level].costEstimator.powerSums[order][1])
                            self.indices[level].costEstimator.powerSums[order][0] = globalIndexPowerSum1
                            self.indices[level].costEstimator.powerSums[order][1] = globalIndexPowerSum2
                            del(globalIndexPowerSum1,globalIndexPowerSum2)
                self.indices[level].costEstimator._sampleCounter = self.indices[level].costEstimator._sampleCounter + self.batchIndices[batch][level].costEstimator._sampleCounter
            # Update model coefficients for cost, bias, variance with new observations
            self.updatePredictors()

            # Update global combined power sums
            for qoi_index in range (self.batchIndices[batch][level].numberOfSamplerOutputs(),self.batchIndices[batch][level].numberOfSamplerOutputs()+self.batchIndices[batch][level].numberOfSamplerCombinedOutputs()):
                if (batch == 0):
                    self.indices[level].qoiEstimator[qoi_index].powerSums = self.batchIndices[batch][level].qoiEstimator[qoi_index].powerSums
                    self.indices[level].qoiEstimator[qoi_index]._sampleCounter = [self.batchIndices[batch][level].qoiEstimator[qoi_index]._sampleCounter]
                else:
                    # asynchronous MC
                    if (self.indices[level].qoiEstimator[qoi_index].order == 1 and self.indices[level].qoiEstimator[qoi_index].indexSetDimension == 0 and self.indices[level].qoiEstimator[qoi_index]):
                        globalIndexPowerSum1,globalIndexPowerSum2 = mcsua.updateGlobalMonteCarloIndexTimePowerSumsOrder2Dimension0_Task(\
                            self.indices[level].qoiEstimator[qoi_index].powerSums[0][0],\
                            self.indices[level].qoiEstimator[qoi_index].powerSums[1][0],\
                            self.batchIndices[batch][level].qoiEstimator[qoi_index].powerSums[0][0],\
                            self.batchIndices[batch][level].qoiEstimator[qoi_index].powerSums[1][0])
                        self.indices[level].qoiEstimator[qoi_index].powerSums[0][0] = globalIndexPowerSum1
                        self.indices[level].qoiEstimator[qoi_index].powerSums[1][0] = globalIndexPowerSum2
                        del(globalIndexPowerSum1,globalIndexPowerSum2)
                    # asynchronous MLMC
                    elif (self.indices[level].qoiEstimator[qoi_index].order == 1 and self.indices[level].qoiEstimator[qoi_index].indexSetDimension == 1 and self.indices[level].qoiEstimator[qoi_index]):
                        globalIndexPowerSumUpper1,globalIndexPowerSumLower1,globalIndexPowerSumUpper2,globalIndexPowerSumLower2 = mcsua.updateGlobalMonteCarloIndexTimePowerSumsOrder2Dimension1_Task(\
                            self.indices[level].qoiEstimator[qoi_index].powerSums[0][0],\
                            self.indices[level].qoiEstimator[qoi_index].powerSums[0][1],\
                            self.indices[level].qoiEstimator[qoi_index].powerSums[1][0],\
                            self.indices[level].qoiEstimator[qoi_index].powerSums[1][1],\
                            self.batchIndices[batch][level].qoiEstimator[qoi_index].powerSums[0][0],\
                            self.batchIndices[batch][level].qoiEstimator[qoi_index].powerSums[0][1],\
                            self.batchIndices[batch][level].qoiEstimator[qoi_index].powerSums[1][0],\
                            self.batchIndices[batch][level].qoiEstimator[qoi_index].powerSums[1][1])
                        # update
                        self.indices[level].qoiEstimator[qoi_index].powerSums[0][0] = globalIndexPowerSumUpper1
                        self.indices[level].qoiEstimator[qoi_index].powerSums[0][1] = globalIndexPowerSumLower1
                        self.indices[level].qoiEstimator[qoi_index].powerSums[1][0] = globalIndexPowerSumUpper2
                        self.indices[level].qoiEstimator[qoi_index].powerSums[1][1] = globalIndexPowerSumLower2
                        del(globalIndexPowerSumUpper1,globalIndexPowerSumLower1,globalIndexPowerSumUpper2,globalIndexPowerSumLower2)
                    else:
                        raise Exception ("Order or dimension not supported. Add the required task into xmc/methodDefs_monteCarloSampler/asynchronousUpdate.py")
                    # append number of samples, since it is a future object and is neeeded only for post processing
                    self.indices[level].qoiEstimator[qoi_index]._sampleCounter.append(self.batchIndices[batch][level].qoiEstimator[qoi_index]._sampleCounter)

    def asynchronousUpdate(self,newHierarchy):
        self.asynchronousInitialize(newHierarchy)
        self.asynchronousLaunchEpoch(newHierarchy)