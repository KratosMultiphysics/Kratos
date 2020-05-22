# XMC imports
import xmc
from xmc.tools import dynamicImport
from xmc.tools import instantiateObject
from xmc.tools import convertObjectToFuture

# Import external libraries
import itertools as it
import numpy as np
import time
import sys
from math import ceil

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
        self.costEstimator = instantiateObject(keywordArgs.get('costEstimator'),**keywordArgs.get('costEstimatorInputDict'))
        # TODO 'qoiEstimator' should be replaced with 'qoiEstimators'.
        # It should also become a list of lists also but this will break
        # things elsewhere in the code and should be done separately later
        self.qoiEstimator = []
        # qoi estimators
        qoi_estimator_module = keywordArgs.get('qoiEstimator')
        qoi_estimator_module_args = keywordArgs.get('qoiEstimatorInputDict')
        for i in range(len(qoi_estimator_module)):
            self.qoiEstimator.append(instantiateObject(qoi_estimator_module[i],**qoi_estimator_module_args[i]))
        # combined estimators
        combined_estimator_module = keywordArgs.get('combinedEstimator')
        combined_estimator_module_args = keywordArgs.get('combinedEstimatorInputDict')
        for i in range(len(combined_estimator_module)):
            self.qoiEstimator.append(instantiateObject(combined_estimator_module[i],**combined_estimator_module_args[i]))

        sampler_input_dict = keywordArgs.get('samplerInputDict')
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
            return self.costEstimator.sampleNumber()

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
        samples = [[] for _ in range (self.numberOfSamplerOutputs()+self.numberOfSamplerCombinedOutputs())]
        time = []
        for _ in range(number_new_samples):
            # Generate a new sample (list of solver outputs for one random event)
            new_sample, new_time = self.newSample()
            # Add each solver output to the correct sublist
            for output_counter in range (self.numberOfSamplerOutputs()+self.numberOfSamplerCombinedOutputs()):
                samples[output_counter].append(new_sample[output_counter])
            time.append([new_time])

        # Pass new samples to relevant estimator
        for output_counter in range(self.numberOfSamplerOutputs()+self.numberOfSamplerCombinedOutputs()):
            # Update estimators for a single solver output
            for mini_batch_counter in range(ceil(len(samples[output_counter])/self.estimatorBatchSize)):
                # Update by estimatorBatchSize
                self.qoiEstimator[output_counter].update(samples[output_counter][mini_batch_counter*self.estimatorBatchSize:(mini_batch_counter+1)*self.estimatorBatchSize])

        # Pass new samples to relevant estimator
        for mini_batch_counter in range(ceil(len(time)/self.estimatorBatchSize)):
            # Update by estimatorBatchSize
            self.costEstimator.update(time[mini_batch_counter*self.estimatorBatchSize:(mini_batch_counter+1)*self.estimatorBatchSize])

    def numberOfSamplerOutputs(self):
        """
        Retrieve the number of outputs from the sampler
        """
        return self.sampler.numberOfOutputs()

    def numberOfSamplerCombinedOutputs(self):
        """
        Retrieve the number of outputs from the sampler
        """
        return self.sampler.numberOfCombinedOutputs()
