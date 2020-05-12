# XMC imports
from xmc.tools import instantiateObject

# Import external libraries
import itertools as it
from math import ceil

import xmc.methodDefs_monteCarloIndex.updateEstimators as mdu

# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

class MonteCarloIndex():
    """
    This is a generic class for an index of an MXMC method. It handles sample
    generation and index-specific estimators.
    """

    # TODO - this constructor should change based on the json file I/O mechanism
    # that will get set up at a future date
    def __init__(self,**keywordArgs):
        # Attributes
        self.indexValue = keywordArgs.get('indexValue')
        self.costEstimator = instantiateObject(keywordArgs.get('costEstimator'),**keywordArgs.get('costEstimatorInputDictionary'))
        # TODO 'qoiEstimator' should be replaced with 'qoiEstimators'.
        # It should also become a list of lists also but this will break
        # things elsewhere in the code and should be done separately later
        self.qoiEstimator = []
        # qoi estimators
        qoi_estimator_module = keywordArgs.get('qoiEstimator')
        qoi_estimator_module_args = keywordArgs.get('qoiEstimatorInputDictionary')
        for i in range(len(qoi_estimator_module)):
            self.qoiEstimator.append(instantiateObject(qoi_estimator_module[i],**qoi_estimator_module_args[i]))
        # combined estimators
        combined_estimator_module = keywordArgs.get('combinedEstimator')
        combined_estimator_module_args = keywordArgs.get('combinedEstimatorInputDictionary')
        for i in range(len(combined_estimator_module)):
            self.qoiEstimator.append(instantiateObject(combined_estimator_module[i],**combined_estimator_module_args[i]))

        sampler_input_dict = keywordArgs.get('samplerInputDictionary')
        sampler_input_dict['solverWrapperIndices'] = self.solverIndices()
        self.sampler = instantiateObject(keywordArgs.get('sampler'),**sampler_input_dict)
        self.samples = keywordArgs.get('samples',None)
        self.areSamplesStored = keywordArgs.get('areSamplesStored',False)
        self.areSamplesRecycled = keywordArgs.get('areSamplesRecycled',True)
        self.estimatorBatchSize = keywordArgs.get('estimatorBatchSize',1)

    def sampleNumber(self):
        if self.areSamplesStored:
            return len(self.samples)
        else:
            return self.qoiEstimator[0].sampleNumber() # necessary for minimizing synchronization points

    # TODO - potentially better name here
    def indexwiseContribution(self, qoiEstimatorCoordinate, qoiEstimatorValueArguements):
        return self.qoiEstimator[qoiEstimatorCoordinate].value(*qoiEstimatorValueArguements)

    def costMoment(self, costEstimatorValueArguments):
        """
        Returns the requested moment estimation of the cost. The input arguments are passed to the cost estimator's relevant method.
        """
        return self.costEstimator.value(*costEstimatorValueArguments)

    def newSample(self):
        """
        Request a new set of correlated samples. Calls its sampler as many times as necessary.
        """
        # Return 2^d QoI vectors that are generated from correlated underlying samples
        list_of_qoi = self.sampler.generate(self.indexValue)
        return list_of_qoi

    def indexSetDimension(self):
        return len(self.indexValue)

    def solverIndices(self):
        """
        Compute list of indices of correlated solvers.
        """
        diff_tuple = []
        hypercube_vertices = 2**self.indexSetDimension()
        for _ in range(self.indexSetDimension()):
            diff_tuple.append([0,1])
        diff_list = list(it.product(*diff_tuple))
        for i in range(len(diff_list)):
            diff_list[i] = list(diff_list[i])

        main_index = [self.indexValue]*hypercube_vertices
        return [[main_index[i][j]-diff_list[i][j] for j in range(self.indexSetDimension())] for i in range(hypercube_vertices)]

    def _addSamples(self,newSamples):
        """
        Add newsamples to stored samples. This is an internal method with no check on whether sample storage is enabled; be careful when using it.
        """
        self.samples.append(newSamples)

    def update(self,newIndexAndSampleNumbers):
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
        if(self.areSamplesRecycled is True):
            number_new_samples = max(5, newIndexAndSampleNumbers[1] - self.sampleNumber())
        else:
            number_new_samples = newIndexAndSampleNumbers[1]
            for i in range(len(self.qoiEstimator)):
                self.qoiEstimator[i].reset()
            self.costEstimator.reset()

        # Generate the required number of correlated samples
        # and estimate cost
        samples = []
        time = []
        for _ in range(number_new_samples):
            # Generate a new sample (list of solver outputs for one random event)
            new_sample, new_time = self.newSample()
            # append to corresponding list
            samples.append(new_sample)
            time.append([new_time])
        # samples = [sample_1, sample_2, .. ]
        # where
        # sample_i = [ [ [QoI_1^l,...,QoI_outputBatchSize^l],[QoI_outputBatchSize+1^l,...,QoI_2*outputBatchSize^l], .. ],
        #            [ [QoI_1^l-1,...,QoI_outputBatchSize^l-1],[QoI_outputBatchSize+1^l-1,...,QoI_2*outputBatchSize^l-1],
        #            .. ] ]

        #Create qoi_groups
        qoi_group=[]
        for mini_qoi_counter in range(ceil(len(self.qoiEstimator)/self.sizeOfSamplerQoiLists())):
            qoi_group.append(self.qoiEstimator[mini_qoi_counter*self.sizeOfSamplerQoiLists():(mini_qoi_counter+1)*self.sizeOfSamplerQoiLists()])

        # Update qoiEstimators groups
        for mini_batch_counter in range(ceil(len(samples)/self.estimatorBatchSize)): # loop samples according to estimatorBatchSize (ebs)
            samples_ebs = samples[mini_batch_counter*self.estimatorBatchSize:(mini_batch_counter+1)*self.estimatorBatchSize]
            for mini_qoi_counter in range(ceil(len(self.qoiEstimator)/self.sizeOfSamplerQoiLists())): # loop qoiEstimators according to outputBatchSize (sql)
                samples_sql = [ [] for _ in range (0,len(samples_ebs))] # obs: len(samples_ebs) = estimatorBatchSize
                count = 0
                for sample_ebs in samples_ebs: # loop over estimatorBatchSize samples
                    for smpls in sample_ebs: # loop over levels of each sample_ebs
                        smpl = smpls[mini_qoi_counter] # corresponding future per level per sample_ebs
                        samples_sql[count].append(smpl)
                    count = count + 1
                mdu.updatePartialQoiEstimators_Task(samples_sql,qoi_group[mini_qoi_counter])
                # delete COMPSs future objects no longer needed
                for smpl in samples_sql:
                    delete_object(*smpl)

        #merge and delete partial updates
        for mini_qoi_counter in range(ceil(len(self.qoiEstimator)/self.sizeOfSamplerQoiLists())):
            mdu.mergeQoiEstimators_Task(qoi_group[mini_qoi_counter], self.qoiEstimator, mini_qoi_counter, self.sizeOfSamplerQoiLists())
            delete_object(qoi_group[mini_qoi_counter])

        # Update costEstimator
        for mini_batch_counter in range(ceil(len(samples)/self.estimatorBatchSize)): # according to estimatorBatchSize (ebs)
            time_ebs = time[mini_batch_counter*self.estimatorBatchSize:(mini_batch_counter+1)*self.estimatorBatchSize]
            mdu.updateCostEstimator_Task(time_ebs,self.costEstimator)

        # delete COMPSs future objects no longer needed
        for t in time:
            delete_object(*t)

    def numberOfSamplerOutputs(self):
        """
        Retrieve the number of outputs from the sampler
        """
        return self.sampler.numberOfOutputs()

    def sizeOfSamplerQoiLists(self):
        """
        Retrieve the number of outputs from the sampler
        """
        return self.sampler.sizeOfQoiLists()
