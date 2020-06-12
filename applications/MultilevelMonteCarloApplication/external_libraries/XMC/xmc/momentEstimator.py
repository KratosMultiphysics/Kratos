# Import xmc classes
import xmc.statisticalEstimator as se
from xmc.tools import dynamicImport

# import python libraries
from math import *

# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

class MomentEstimator(se.StatisticalEstimator):
    """
    This class estimates raw and central moments up to a given order (including
    expectations, obviously). These estimations are computed using h-statistics,
    so we store (and update) the power sums from which they are computed.
    """

    def __init__(self,**keywordArgs):
        # Parent constructor
        super().__init__(**keywordArgs)

        # Attributes
        # maximum order to which the moments are computed
        self.order = keywordArgs.get('order')
        # TODO indexSetDimension is only used to initialise self.powerSums
        # and define the default value of self._updatedPowerSums.
        # Should it be a required input argument?
        indexSetDimension = keywordArgs.get('indexSetDimension')
        self.indexSetDimension = indexSetDimension
        # list of power sums required to compute the moments
        # TODO - This is a temporary fix until when sample_moments.py is interfaced
        if indexSetDimension==0:
            self.powerSums = [[None] for _ in range(self.powerSumsOrder())]
        elif indexSetDimension==1:
            self.powerSums = [[None for _ in range(i+2)]
                              for i in range (self.powerSumsOrder())]
        else:
            self.powerSums = None

        # Methods
        basePath = 'xmc.methodDefs_momentEstimator.'
        self._centralMomentComputer = dynamicImport(
            keywordArgs.get('centralMomentComputer', basePath+
                            'computeCentralMoments.centralMomentWrapper'))
        self._centralMomentErrorComputer = dynamicImport(
            keywordArgs.get('centralMomentErrorComputer', basePath+
                            'computeErrorEstimation.centralMomentErrorWrapper'))
        self._rawMomentComputer = dynamicImport(keywordArgs.get('rawMomentComputer',None))
        self._rawMomentErrorComputer = dynamicImport(
            keywordArgs.get('rawMomentErrorComputer',None))
        # updatedPowerSums method (conditional default value)
        defaultUpdater = (basePath + 'updatePowerSums.updatePowerSums'
                          'Order{o}Dimension{d}').format(
                              o = self.powerSumsOrder(), d = self.indexSetDimension)
        self._updatedPowerSums = dynamicImport(
            keywordArgs.get('updatedPowerSums', defaultUpdater))

    def update(self,samples):
        """
        function updating the power sums given new samples
        input:  self:    an instance of the class
        samples: list containing new samples
        """
        
        dimension = log2(len(samples[0])) # len(power_sum[0])=2**dimension
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


    def value(self,order,isCentral=None,isErrorEstimationRequested=False):
        """
        function returning a specific order moment, its error estimate or both
        Inputs:
        self: an instance of the class
        order: order of the moment
        isCentral: setting if the working moment is central or raw
        computeErrorEstimation: setting if computing or not the errorEstimation
        """
        
        if isCentral is None:
            isCentral = order>1
        # TODO the [] around the returns in the following is to ensure consistency with
        # the inputs expected by the rest of the program. This is becuase the quantity
        # returned can be either a scalar estimation value or a list of values of function
        # evalutaions at certain points in a domain (eg. CDF)
        if isCentral:
            return self._computeCentralMoment(order,isErrorEstimationRequested)
        else:
            return self._computeRawMoment(order,isErrorEstimationRequested)

    def errorEstimation(self,order,isCentral):
        """
        function returning the error estimation of working moment from stored power sums
        input:  self: an instance of the class
            order: order of the moment
        """
        
        #TODO not normally used any more. Keep ?
        #TODO - before inferring dimension from powerSums, need to ensure that
        # update method has been run at least a certain number of times required
        # to compute the moment of requested order
        dimension = log2(len(self.powerSums[0])) # len(power_sum[0])=2**dimension
        # Let's unpack the nested list self.powerSums
        # First all elements in self.powerSums[0], then those in self.powerSums[1], etc.
        # TODO is it the order in which we want this list?
        power_sums = [item for sublist in self.powerSums for item in sublist]
        if isCentral:
            return self._centralMomentErrorComputer(
                dimension,order,*[*power_sums,self.sampleNumber()])
        else:
            return self._rawMomentErrorComputer(
                dimension,order,*[*power_sums,self.sampleNumber()])

    def _computeCentralMoment(self,order,isErrorEstimationRequested):
        """
        function returning the central moment of working moment from stored power sums
        input:  self: an instance of the class
            order: order of the moment
        """
        
        dimension = log2(len(self.powerSums[0])) # len(power_sum[0])=2**dimension
        # Let's unpack the nested list self.powerSums
        # First all elements in self.powerSums[0], then those in self.powerSums[1], etc.
        # TODO is it the order in which we want this list?
        if isErrorEstimationRequested:
            power_sums = [item for sublist in self.powerSums[:2*order] for item in sublist]
            output = self._centralMomentErrorComputer(
                dimension,order,*[*power_sums,self.sampleNumber()])
        else:
            power_sums = [item for sublist in self.powerSums[:order] for item in sublist]
            output =  self._centralMomentComputer(
                dimension,order,*[*power_sums,self.sampleNumber()])
        return output

    def _computeRawMoment(self,order,isErrorEstimationRequested):
        """
        Compute raw moment of requested order from stored power sums.
        """
        
        dimension = log2(len(self.powerSums[0])) # len(power_sum[0])=2**dimension
        # Let's unpack the nested list self.powerSums
        # First all elements in self.powerSums[0], then those in self.powerSums[1], etc.
        # TODO is it the order in which we want this list?
        power_sums = [item for sublist in self.powerSums for item in sublist]
        output = self._rawMomentComputer(dimension,order,*power_sums,self.sampleNumber())
        if isErrorEstimationRequested:
            output = [output, self._rawMomentErrorComputer(
                dimension,order,*power_sums,self.sampleNumber())]
        return output

    def powerSumsOrder(self):
        """
        Returns the maximum order to which the powerSums are computed to compute estimation and error of moments uo to order self.order.
        """
        
        return 2*self.order

    def reset(self):
        """
        Reset the value of power sums to zero. Reset the number of samples
        """
        
        noneList = [[None]*len(powerSum) for powerSum in self.powerSums]
        self.powerSums = noneList
        self._sampleCounter = 0


class CombinedMomentEstimator(MomentEstimator):
    """
    This class estimates raw and central moments up to a given order (including
    expectations, obviously). These estimations are computed using h-statistics,
    so we store (and update) the power sums from which they are computed.
    The computed statistics are the composition of sampling and time statistical processes.
    """

    def __init__(self,**keywordArgs):
        # Parent constructor
        super(CombinedMomentEstimator,self).__init__(**keywordArgs)

        # list of power sums to combined compute moments
        # for these power sums we store upper and lower level
        # for example, for index = 1 and order = 1, we will have S10, S01, S20, S02
        if (self.indexSetDimension == 0):
            self.powerSums = [[None] for _ in range (self.powerSumsOrder())]
        elif (self.indexSetDimension == 1):
            self.powerSums = [[None for _ in range(0,2)]
                              for _ in range (self.powerSumsOrder())]
        else:
            self.powerSums = None

        # Methods
        basePath = 'xmc.methodDefs_momentEstimator.'
        self._centralMomentComputer = dynamicImport(
            keywordArgs.get('centralMomentComputer', basePath+
            'computeCentralMoments.centralCombinedMomentWrapper'))
        self._centralMomentErrorComputer = dynamicImport(
            keywordArgs.get('centralMomentErrorComputer', basePath+
            'computeErrorEstimation.centralCombinedMomentErrorWrapper'))
        self._rawMomentComputer = dynamicImport(
            keywordArgs.get('rawMomentComputer',None))
        self._rawMomentErrorComputer = dynamicImport(
            keywordArgs.get('rawMomentErrorComputer',None))
        # updatedPowerSums method (conditional default value)
        defaultUpdater = (basePath + 'updateCombinedPowerSums.updatePowerSums'
                          'Order{o}Dimension{d}').format(
                              o=self.powerSumsOrder(), d=self.indexSetDimension)
        self._updatedPowerSums = dynamicImport(
            keywordArgs.get('updatedTimePowerSums', defaultUpdater))

    def update(self,samples):
        """
        function updating the power sums given new samples
        input:  self:    an instance of the class
        samples: list containing new samples
        """
        # Let's unpack the nested list self.powerSums
        # First all elements in self.powerSums[0], then those in self.powerSums[1], etc.
        # Idem for samples
        # TODO is it the order in which we want these lists?
        power_sums = [item for sublist in self.powerSums for item in sublist]
        sample_counter_power_sums = self._updatedPowerSums(self._sampleCounter,samples,*power_sums)
        sample_counter = sample_counter_power_sums[0]
        power_sums = list(sample_counter_power_sums[1:])
        for i in range(len(self.powerSums)):
            for j in range(len(self.powerSums[i])):
                self.powerSums[i][j] = power_sums.pop(0)
        self._sampleCounter = sample_counter
