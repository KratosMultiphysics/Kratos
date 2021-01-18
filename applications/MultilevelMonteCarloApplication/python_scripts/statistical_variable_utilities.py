from __future__ import absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import packages
import numpy as np

# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3


"""
auxiliary function of UpdateOnePassCentralMoments of the StatisticalVariable class
input:  sample: new value that will update the statistics
        old_mean             : old mean
        old_central_moment_1 : old first central moment
        compute_M1           : boolean setting if computation is needed
        old_central_moment_2 : old second central moment
        compute_M2           : boolean setting if computation is needed
        old_central_moment_3 : old third central moment
        compute_M3           : boolean setting if computation is needed
        old_central_moment_1 : old fourth central moment
        compute_M4           : boolean settings if computation is needed
        nsamples             : old number of samples computed, starts from 1
output: new_mean             : updated mean
        new_sample_variance  : updated sample variance
        new_central_moment_1 : updated central_moment_1
        new_central_moment_2 : updated central_moment_2
        new_central_moment_3 : updated central_moment_3
        new_central_moment_4 : updated central_moment_4
        nsamples             : updated number of samples
"""
@constraint(ComputingUnits="${computing_units_auxiliar_utilities}")
@ExaquteTask(returns=7,priority=True)
def UpdateOnePassCentralMomentsAux_Task(sample,old_mean,old_central_moment_1,compute_M1,old_central_moment_2,compute_M2,old_central_moment_3,compute_M3,old_central_moment_4,compute_M4,nsamples):
    old_M1 = old_central_moment_1 * nsamples
    old_M2 = old_central_moment_2 * nsamples
    old_M3 = old_central_moment_3 * nsamples
    old_M4 = old_central_moment_4 * nsamples
    nsamples = nsamples + 1
    if nsamples == 1:
        new_mean = sample
        new_M1 = 0.0
        new_M2 = 0.0
        new_sample_variance = 0.0
        new_M3 = 0.0
        new_M4 = 0.0
    else:
        delta = np.subtract(sample,old_mean)
        new_mean = old_mean + np.divide(delta,nsamples)
        if (compute_M1):
            new_M1 = old_M1 # we are not updating, first central moment = 0.0
        else:
            new_M1 = old_M1 # we are not updating, first central moment = 0.0
        if (compute_M2):
            new_M2 = old_M2 + delta*np.subtract(sample,new_mean)
        else:
            raise Exception ("Not computing StatisticalVariable.central_moment_2, set StatisticalVariable.central_moment_2_to_compute to True")
        new_sample_variance = np.divide(new_M2,np.subtract(nsamples,1))
        if (compute_M3):
            new_M3 = old_M3 - 3.0*old_M2*np.divide(delta,nsamples) + np.divide(np.multiply((nsamples-1)*(nsamples-2),(delta**3)),(nsamples**2))
        else:
            new_M3 = old_M3 # we are not updating
        if (compute_M4):
            new_M4 = old_M4 - 4.0*old_M3*np.divide(delta,nsamples) + 6.0*old_M2*np.divide(delta,nsamples)**2 + np.multiply((nsamples-1)*(nsamples**2-3*nsamples+3),np.divide(delta**4,nsamples**3))
        else:
            new_M4 = old_M4 # we are not updating
    new_central_moment_1 = new_M1 / nsamples
    new_central_moment_2 = new_M2 / nsamples
    new_central_moment_3 = new_M3 / nsamples
    new_central_moment_4 = new_M4 / nsamples
    return new_mean,new_sample_variance,new_central_moment_1,new_central_moment_2,new_central_moment_3,new_central_moment_4,nsamples


"""
auxiliary function of UpdateOnePassPowerSums of the StatisticalVariable class
input:  sample   : new value that will update the statistics
        old_S1   : old first power sum
        old_S2   : old second power sum
        old_S3   : old third power sum
        old_S4   : old fourth power sum
        nsamples : number of samples, it has already been updated in UpdateOnePassCentralMomentsAux_Task
output: new_S1 : updated first power sum
        new_s2 : updated second power sum
        new_S3 : updated third power sum
        new_S4 : updated fourth power sum
"""
@constraint(ComputingUnits="${computing_units_auxiliar_utilities}")
@ExaquteTask(returns=5,priority=True)
def UpdateOnePassPowerSumsAux_Task(sample,old_S1,old_S2,old_S3,old_S4,nsamples):
    nsamples = nsamples + 1
    if nsamples == 1:
        new_S1 = sample
        new_S2 = sample**2
        new_S3 = sample**3
        new_S4 = sample**4
    else:
        new_S1 = old_S1 + sample
        new_S2 = old_S2 + sample**2
        new_S3 = old_S3 + sample**3
        new_S4 = old_S4 + sample**4
    return new_S1,new_S2,new_S3,new_S4,nsamples


"""
auxiliary function of UpdateGlobalPowerSums of the StatisticalVariable class
input:  old_S1                   : old first power sum
        old_S2                   : old second power sum
        old_S3                   : old third power sum
        old_S4                   : old fourth power sum
        number_samples_level     : number of samples, it has already been updated in UpdateOnePassCentralMomentsAux_Task
        add_S1                   : power sum order one to add
        add_S2                   : power sum order two to add
        add_S3                   : power sum order three to add
        add_S4                   : power sum order four to add
        add_number_samples_level : number of samples to add
output: new_S1               : updated first power sum
        new_s2               : updated second power sum
        new_S3               : updated third power sum
        new_S4               : updated fourth power sum
        number_samples_level : number of samples of current level
"""
@constraint(ComputingUnits="${computing_units_auxiliar_utilities}")
@ExaquteTask(returns=5,priority=True)
def UpdateGlobalPowerSumsAux_Task(old_S1,old_S2,old_S3,old_S4,number_samples_level,add_S1,add_S2,add_S3,add_S4,add_number_samples_level):
    new_S1 = old_S1 + add_S1
    new_S2 = old_S2 + add_S2
    new_S3 = old_S3 + add_S3
    new_S4 = old_S4 + add_S4
    number_samples_level = number_samples_level + add_number_samples_level
    return new_S1,new_S2,new_S3,new_S4,number_samples_level


"""
function unfolding values from a list, needed by PyCOMPSs for list of lists
input:  sample  : the list of lists
output: sample[*] : list position * of the list of lists
"""
@constraint(ComputingUnits="${computing_units_auxiliar_utilities}")
@ExaquteTask(returns=4, priority=True)
def UnfoldValuesAux_Task(sample):
    return sample[0], sample[1], sample[2], sample[3]


"""
auxiliary function of UpdateBatchesPassPowerSums
input:  samples : list of samples
output: return the sum, done in mini_batch_size batches, of the samples components
"""
@constraint(ComputingUnits="${computing_units_auxiliar_utilities}")
@ExaquteTask(returns=1,priority=True)
def UpdateBatchesPassPowerSumsAux_Task(*samples):
    samples_list = np.array(list(samples))
    return np.sum(samples_list, axis = 0)
    """
    if nsamples == 0:
        new_S1 = samples[0]
        new_S2 = samples[0]**2
        new_S3 = samples[0]**3
        new_S4 = samples[0]**4
        old_S1 = new_S1
        old_S2 = new_S2
        old_S3 = new_S3
        old_S4 = new_S4
        nsamples = 1
        samples=samples[1:]
    for sample in samples:
        nsamples = nsamples + 1
        new_S1 = old_S1 + sample
        new_S2 = old_S2 + sample**2
        new_S3 = old_S3 + sample**3
        new_S4 = old_S4 + sample**4
        old_S1 = new_S1
        old_S2 = new_S2
        old_S3 = new_S3
        old_S4 = new_S4
    return new_S1,new_S2,new_S3,new_S4,nsamples
    """

"""
auxiliary function of UpdateHStatistics of the StatisticalVariable class
input:  S1_level             : first power sum at defined level
        S2_level             : second power sum at defined level
        S3_level             : third power sum at defined level
        S4_level             : fourth power sum at defined level
        number_samples_level : number of samples (already update) for defined level
output: h1_level : first h statistics for defined level
        h2_level : second h statistics for defined level
        h3_level : third h statistics for defined level
        h4_level : fourth h statistics for defined level
"""
@constraint(ComputingUnits="${computing_units_auxiliar_utilities}")
@ExaquteTask(returns=4,priority=True)
def ComputeHStatisticsAux_Task(S1_level,S2_level,S3_level,S4_level,number_samples_level):
    h1_level = S1_level / number_samples_level
    h2_level = (number_samples_level*S2_level-S1_level**2) / ((number_samples_level-1)*number_samples_level)
    h3_level = (number_samples_level**2*S3_level-3*number_samples_level*S2_level*S1_level+2*S1_level**3) / \
        ((number_samples_level-2)*(number_samples_level-1)*number_samples_level)
    h4_level = ((-4*number_samples_level**2+8*number_samples_level-12)*S3_level*S1_level+ \
        (number_samples_level**3-2*number_samples_level**2+3*number_samples_level)*S4_level+ \
        6*number_samples_level*S2_level*S1_level**2+(9-6*number_samples_level)*S2_level**2-3*S1_level**4) / \
        ((number_samples_level-3)*(number_samples_level-2)*(number_samples_level-1)*number_samples_level)
    return h1_level,h2_level,h3_level,h4_level


"""
auxiliary function of ComputeSkewnessKurtosis of the StatisticalVariable class
input:  h2_level : second h statistics for defined level
        h3_level : third h statistics for defined level
        h4_level : fourth h statistics for defined level
output: skewness_level : skewness for defined level
        kurtosis_level : kurtosis for defined level
"""
@constraint(ComputingUnits="${computing_units_auxiliar_utilities}")
@ExaquteTask(returns=2,priority=True)
def ComputeSkewnessKurtosisAux_Task(h2_level,h3_level,h4_level):
    skewness_level = h3_level / (np.sqrt(h2_level**3))
    kurtosis_level = h4_level / (h2_level**2)
    return skewness_level,kurtosis_level


"""
auxiliary function of ComputeSampleCentralMomentsFromScratch of the StatisticalVariable class
input:  sample: new value that will update the statistics
        number_samples_level                              : number of samples for defined level
        central_moment_from_scratch_1_to_compute          : boolean setting if computation is needed
        central_moment_from_scratch_2_to_compute          : boolean setting if computation is needed
        central_moment_from_scratch_3_to_compute          : boolean setting if computation is needed
        central_moment_from_scratch_3_absolute_to_compute : boolean setting if computation is needed
        central_moment_from_scratch_4_to_compute          : boolean setting if computation is needed
        central_moment_from_scratch_1                     : old first central moment
        central_moment_from_scratch_2                     : old second central moment
        central_moment_from_scratch_3                     : old third central moment
        central_moment_from_scratch_3_absolute            : old third central moment absolute value
        central_moment_from_scratch_4                     : old fourth central moment
output: central_moment_from_scratch_1          : updated first central moment
        central_moment_from_scratch_2          : updated second central moment
        central_moment_from_scratch_3          : updated third central moment
        central_moment_from_scratch_3_absolute : updated third central moment absolute value
        central_moment_from_scratch_4          : update fourth central moment
"""
@constraint(ComputingUnits="${computing_units_auxiliar_utilities}")
@ExaquteTask(returns=5,priority=True)
def ComputeSampleCentralMomentsFromScratchAux_Task(number_samples_level,central_moment_from_scratch_1_to_compute,central_moment_from_scratch_2_to_compute, \
    central_moment_from_scratch_3_to_compute,central_moment_from_scratch_3_absolute_to_compute,central_moment_from_scratch_4_to_compute, \
    central_moment_from_scratch_1,central_moment_from_scratch_2,central_moment_from_scratch_3,central_moment_from_scratch_3_absolute,central_moment_from_scratch_4, \
    samples):
    # generate a single list from a list of lists
    samples = [item for sublist in samples for item in sublist]
    # compute the mean
    auxiliary_mean = 0.0
    for sample in samples:
        auxiliary_mean = auxiliary_mean + sample
    curr_mean = auxiliary_mean / number_samples_level
    for sample in samples:
        if (central_moment_from_scratch_1_to_compute):
            central_moment_from_scratch_1 = central_moment_from_scratch_1 + ((sample - curr_mean)**1) / number_samples_level
        if (central_moment_from_scratch_2_to_compute):
            central_moment_from_scratch_2 = central_moment_from_scratch_2 + ((sample - curr_mean)**2) / number_samples_level
        if (central_moment_from_scratch_3_to_compute):
            central_moment_from_scratch_3 = central_moment_from_scratch_3 + ((sample - curr_mean)**3) / number_samples_level
        if (central_moment_from_scratch_3_absolute_to_compute):
            central_moment_from_scratch_3_absolute = central_moment_from_scratch_3_absolute + (np.abs(sample - curr_mean)**3) / number_samples_level
        if (central_moment_from_scratch_4_to_compute):
            central_moment_from_scratch_4 = central_moment_from_scratch_4 + ((sample - curr_mean)**4) / number_samples_level
    return central_moment_from_scratch_1,central_moment_from_scratch_2,central_moment_from_scratch_3,central_moment_from_scratch_3_absolute,central_moment_from_scratch_4


class StatisticalVariable(object):
    """The base class for statistical variables"""
    def __init__(self):
        """constructor of the class
        Keyword arguments:
        self : an instance of a class
        """

        # values of the variable, organized per level
        self.values = []
        # mean of the variable per each level
        self.raw_moment_1 = []
        # sample variance of the variable per each level
        self.unbiased_central_moment_2 = []
        # moments of the variable per each level M_p  = n * mu_p
        #                                        mu_p = p-th central moment
        #                                        n    = number of values
        self.central_moment_1 = []
        self.central_moment_2 = []
        self.central_moment_3 = []
        self.central_moment_4 = []
        # set which central moments will be computed (moment_2 is mandatory to be computed because it is exploited in the mean evaluation)
        self.central_moment_1_to_compute = True
        self.central_moment_2_to_compute = True
        self.central_moment_3_to_compute = True
        self.central_moment_4_to_compute = True
        # bias error of the variable
        self.bias_error = None
        # statistical error of the variable
        self.statistical_error = None
        # type of variable: scalar or field
        self.type = None
        # number of samples of the variable
        self.number_samples = None
        self.batches_number_samples = []
        # global power sums
        # S_p = \sum_{i=1}^{n} Q(sample_i)**p, organized per level
        self.power_sum_1 = []
        self.power_sum_2 = []
        self.power_sum_3 = []
        self.power_sum_4 = []
        # power sums batches
        self.power_sum_batches_1 = []
        self.power_sum_batches_2 = []
        self.power_sum_batches_3 = []
        self.power_sum_batches_4 = []
        # sample central moments \mu_p = \sum_{i=1}^{n} (Q(sample_i)-mean_n)**p / n, organized per level
        self.central_moment_from_scratch_1 = []
        self.central_moment_from_scratch_2 = []
        self.central_moment_from_scratch_3 = []
        self.central_moment_from_scratch_3_absolute = [] # \mu_p = \sum_{i=1}^{n} abs((Q(sample_i)-mean_n)**p) / n
        self.central_moment_from_scratch_4 = []
        self.central_moment_from_scratch_1_to_compute = False
        self.central_moment_from_scratch_2_to_compute = False
        self.central_moment_from_scratch_3_to_compute = False
        self.central_moment_from_scratch_3_absolute_to_compute = False
        self.central_moment_from_scratch_4_to_compute = False
        # h-statistics h_p, the unbiased central moment estimator with minimal variance, organized per level
        self.h_statistics_1 = []
        self.h_statistics_2 = []
        self.h_statistics_3 = []
        self.h_statistics_4 = []
        self.h_statistics_computed = False
        # skewness of the variable per each level
        self.skewness = []
        # kurtosis of the variable per each level
        self.kurtosis = []
        # convergence criteria of the algorithm
        self.convergence_criteria = None

    """
    function initializing variables of the Statistical Variable class in lists given number of levels
    input:  self                   : an instance of the class
            number_levels          : number of levels considered
            number_initial_batches : number of batches of iteration zero
    """
    def InitializeLists(self,number_levels,number_initial_batches):
        self.number_samples = [0 for _ in range (number_levels)]
        self.values = [[[] for _ in range (number_levels)] for _ in range (number_initial_batches)]
        self.raw_moment_1 = [[] for _ in range (number_levels)]
        self.central_moment_1 = [[] for _ in range (number_levels)]
        self.central_moment_2 = [[] for _ in range (number_levels)]
        self.central_moment_3 = [[] for _ in range (number_levels)]
        self.central_moment_4 = [[] for _ in range (number_levels)]
        self.unbiased_central_moment_2 = [[] for _ in range (number_levels)]
        self.power_sum_1 = [0 for _ in range (number_levels)]
        self.power_sum_2 = [0 for _ in range (number_levels)]
        self.power_sum_3 = [0 for _ in range (number_levels)]
        self.power_sum_4 = [0 for _ in range (number_levels)]
        self.power_sum_batches_1 = [[[] for _ in range (number_levels)] for _ in range (number_initial_batches)]
        self.power_sum_batches_2 = [[[] for _ in range (number_levels)] for _ in range (number_initial_batches)]
        self.power_sum_batches_3 = [[[] for _ in range (number_levels)] for _ in range (number_initial_batches)]
        self.power_sum_batches_4 = [[[] for _ in range (number_levels)] for _ in range (number_initial_batches)]
        self.h_statistics_1 = [[] for _ in range (number_levels)]
        self.h_statistics_2 = [[] for _ in range (number_levels)]
        self.h_statistics_3 = [[] for _ in range (number_levels)]
        self.h_statistics_4 = [[] for _ in range (number_levels)]
        self.skewness = [[] for _ in range (number_levels)]
        self.kurtosis = [[] for _ in range (number_levels)]
        self.central_moment_from_scratch_1 = [[] for _ in range (number_levels)]
        self.central_moment_from_scratch_2 = [[] for _ in range (number_levels)]
        self.central_moment_from_scratch_3 = [[] for _ in range (number_levels)]
        self.central_moment_from_scratch_3_absolute = [[] for _ in range (number_levels)]
        self.central_moment_from_scratch_4 = [[] for _ in range (number_levels)]
        self.batches_number_samples = [[0 for _ in range (number_levels)] for _ in range (number_initial_batches)]

    """
    function updating statistic moments and number of samples
    input:  self     : an instance of the class
            level    : defined level
            i_sample : defined sample in level
    """
    def UpdateOnePassCentralMoments(self,level,i_sample):
        number_samples_level = self.number_samples[level]
        sample = self.values[level][i_sample]
        old_mean = self.raw_moment_1[level]
        # old_M1 = self.central_moment_1[level] * number_samples_level
        old_central_moment_1 = self.central_moment_1[level]
        compute_M1 = self.central_moment_1_to_compute
        # old_M2 = self.central_moment_2[level] * number_samples_level
        old_central_moment_2 = self.central_moment_2[level]
        compute_M2 = self.central_moment_2_to_compute
        # old_M3 = self.central_moment_3[level] * number_samples_level
        old_central_moment_3 = self.central_moment_3[level]
        compute_M3 = self.central_moment_3_to_compute
        # old_M4 = self.central_moment_4[level] * number_samples_level
        old_central_moment_4 = self.central_moment_4[level]
        compute_M4 = self.central_moment_4_to_compute
        new_mean,new_sample_variance,new_central_moment_1,new_central_moment_2,new_central_moment_3,new_central_moment_4,number_samples_level \
            = UpdateOnePassCentralMomentsAux_Task(sample,old_mean,old_central_moment_1,compute_M1,old_central_moment_2,compute_M2,old_central_moment_3,compute_M3,old_central_moment_4,compute_M4,number_samples_level)
        self.raw_moment_1[level] = new_mean
        self.unbiased_central_moment_2[level] = new_sample_variance
        self.central_moment_1[level] = new_central_moment_1
        self.central_moment_2[level] = new_central_moment_2
        self.central_moment_3[level] = new_central_moment_3
        self.central_moment_4[level] = new_central_moment_4
        self.number_samples[level] = number_samples_level

    """
    function updating the power sums S_p
    input:  self     : an instance of the class
            level    : defined level
            i_sample : defined sample in level
    """
    def UpdateOnePassPowerSums(self,level,i_sample):
        sample = self.values[level][i_sample]
        old_S1 = self.power_sum_1[level]
        old_S2 = self.power_sum_2[level]
        old_S3 = self.power_sum_3[level]
        old_S4 = self.power_sum_4[level]
        number_samples_level = self.number_samples[level]
        new_S1,new_S2,new_S3,new_S4,number_samples_level = UpdateOnePassPowerSumsAux_Task(sample,old_S1,old_S2,old_S3,old_S4,number_samples_level)
        self.power_sum_1[level] = new_S1
        self.power_sum_2[level] = new_S2
        self.power_sum_3[level] = new_S3
        self.power_sum_4[level] = new_S4
        self.number_samples[level] = number_samples_level

    """
    function updating the global power sums
    input:  self          : an instance of the class
            level         : current level
            batch_counter : current batch
    """
    def UpdateGlobalPowerSums(self,level,batch_counter):
        old_S1 = self.power_sum_1[level]
        old_S2 = self.power_sum_2[level]
        old_S3 = self.power_sum_3[level]
        old_S4 = self.power_sum_4[level]
        number_samples_level = self.number_samples[level]

        add_S1 = self.power_sum_batches_1[batch_counter][level]
        add_S2 = self.power_sum_batches_2[batch_counter][level]
        add_S3 = self.power_sum_batches_3[batch_counter][level]
        add_S4 = self.power_sum_batches_4[batch_counter][level]
        add_number_samples_level = self.batches_number_samples[batch_counter][level]

        new_S1,new_S2,new_S3,new_S4,number_samples_level = UpdateGlobalPowerSumsAux_Task(old_S1,old_S2,old_S3,old_S4,number_samples_level,add_S1,add_S2,add_S3,add_S4,add_number_samples_level)
        self.power_sum_1[level] = new_S1
        self.power_sum_2[level] = new_S2
        self.power_sum_3[level] = new_S3
        self.power_sum_4[level] = new_S4
        self.number_samples[level] = number_samples_level

    """
    function updating the in-batch power sums
    input:  self          : an instance of the class
            level         : current level
            batch_counter : current batch
            mini_batch    : size such that we update the in-batch power sums with mini_batch samples
    """
    def UpdateBatchesPassPowerSum(self,level,batch_counter,mini_batch=50):
        samples = self.values[batch_counter][level]
        #for mini_batch in range (0,len(samples)):
        while len(samples) > 1:
            mini_batches_samples = samples[:mini_batch]
            samples = samples[mini_batch:]
            new_power_sums = UpdateBatchesPassPowerSumsAux_Task(*mini_batches_samples)
            samples.append(new_power_sums)
        new_S1, new_S2, new_S3, new_S4 = UnfoldValuesAux_Task(samples[0])
        self.power_sum_batches_1[batch_counter][level] = new_S1
        self.power_sum_batches_2[batch_counter][level] = new_S2
        self.power_sum_batches_3[batch_counter][level] = new_S3
        self.power_sum_batches_4[batch_counter][level] = new_S4


    """
    function computing the h statistics h_p from the power sums
    input:  self  : an instance of the class
            level : defined level
    """
    def ComputeHStatistics(self,level):
        number_samples_level = self.number_samples[level]
        S1_level = self.power_sum_1[level]
        S2_level = self.power_sum_2[level]
        S3_level = self.power_sum_3[level]
        S4_level = self.power_sum_4[level]
        self.h_statistics_computed = True
        h1_level,h2_level,h3_level,h4_level = ComputeHStatisticsAux_Task(S1_level,S2_level,S3_level,S4_level,number_samples_level)
        self.h_statistics_1[level] = h1_level
        self.h_statistics_2[level] = h2_level
        self.h_statistics_3[level] = h3_level
        self.h_statistics_4[level] = h4_level

    """
    function computing from scratch the central moments and the absolute third central moment
    input:  self  : an instance of the class
            level : defined level
    """
    def ComputeSampleCentralMomentsFromScratch(self,level,number_samples_level):
        # initialize central moments
        central_moment_from_scratch_1 = 0.0
        central_moment_from_scratch_2 = 0.0
        central_moment_from_scratch_3 = 0.0
        central_moment_from_scratch_3_absolute = 0.0
        central_moment_from_scratch_4 = 0.0
        central_moment_from_scratch_1_to_compute = self.central_moment_from_scratch_1_to_compute
        central_moment_from_scratch_2_to_compute = self.central_moment_from_scratch_2_to_compute
        central_moment_from_scratch_3_to_compute = self.central_moment_from_scratch_3_to_compute
        central_moment_from_scratch_3_absolute_to_compute = self.central_moment_from_scratch_3_absolute_to_compute
        central_moment_from_scratch_4_to_compute = self.central_moment_from_scratch_4_to_compute

        samples = []
        for batch in range (len(self.values)):
            for mini_batch_samples in self.values[batch][level]:
                samples.append(mini_batch_samples)

        central_moment_from_scratch_1,central_moment_from_scratch_2,central_moment_from_scratch_3,central_moment_from_scratch_3_absolute,central_moment_from_scratch_4 = \
            ComputeSampleCentralMomentsFromScratchAux_Task(number_samples_level,central_moment_from_scratch_1_to_compute, \
            central_moment_from_scratch_2_to_compute,central_moment_from_scratch_3_to_compute,central_moment_from_scratch_3_absolute_to_compute,central_moment_from_scratch_4_to_compute, \
            central_moment_from_scratch_1,central_moment_from_scratch_2,central_moment_from_scratch_3,central_moment_from_scratch_3_absolute,central_moment_from_scratch_4, samples)
        self.central_moment_from_scratch_1[level] = central_moment_from_scratch_1
        self.central_moment_from_scratch_2[level] = central_moment_from_scratch_2
        self.central_moment_from_scratch_3[level] = central_moment_from_scratch_3
        self.central_moment_from_scratch_3_absolute[level] = central_moment_from_scratch_3_absolute
        self.central_moment_from_scratch_4[level] = central_moment_from_scratch_4

    """
    function computing the skewness and the kurtosis from the h statistics
    skewness = \mu_3 / \sqrt(\mu_2^3)
    kurtosis = \mu_4 / \mu_2^2
    input:  self  : an instance of the class
            level : defined level
    """
    def ComputeSkewnessKurtosis(self,level):
        if (self.h_statistics_computed):
            h2_level = self.h_statistics_2[level]
            h3_level = self.h_statistics_3[level]
            h4_level = self.h_statistics_4[level]
            skewness_level,kurtosis_level =ComputeSkewnessKurtosisAux_Task(h2_level,h3_level,h4_level)
            self.skewness[level] = skewness_level
            self.kurtosis[level] = kurtosis_level
