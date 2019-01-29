from __future__ import absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import packages
import numpy as np

# Import exaqute
from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with pycompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
# from exaqute.ExaquteTaskLocal import *      # to execute with python3
'''
get_value_from_remote is the equivalent of compss_wait_on: a synchronization point
in future, when everything is integrated with the it4i team, importing exaqute.ExaquteTaskHyperLoom you can launch your code with their scheduler instead of BSC
'''

'''
This utility contains the classes to perform the Continuation Multilevel Monte Carlo (CMLMC) and the Monte Carlo (MC) algorithms

References:
P. Pébay, T. B. Terriberry, H. Kolla, J. Bennett; Stable, Scalable Formulas for Parallel and Online Computation of Higher-order Multivariate Central Moments with Arbitrary Weights; Comput. Stat.(2016) 31:1305-1325
M. Pisaroni, S. Krumscheid, F. Nobile; Quantifying uncertain system outputs via the multilevel Monte Carlo method - Part I: Central moment estimation; available at MATHICSE technical report no. 23.2017
'''


'''
auxiliary function of UpdateOnePassMomentsVariance of the StatisticalVariable class
'''
@ExaquteTask(returns=6)
def UpdateOnePassMomentsVarianceAux_Task(sample,old_mean,old_M2,old_M3,compute_M3,old_M4,compute_M4,nsamples):
    nsamples = nsamples + 1
    if nsamples == 1:
        new_mean = sample
        new_M2 = np.zeros(np.size(sample))
        new_M2 = np.asscalar(new_M2) # do so to have a list of scalars, and not a list of arrays of one element
        new_sample_variance = np.zeros(np.size(sample))
        new_sample_variance = np.asscalar(new_sample_variance) # do so to have a list of scalars, and not a list of arrays of one element
        new_M3 = np.zeros(np.size(sample))
        new_M3 = np.asscalar(new_M3) # do so to have a list of scalars, and not a list of arrays of one element
        new_M4 = np.zeros(np.size(sample))
        new_M4 = np.asscalar(new_M4) # do so to have a list of scalars, and not a list of arrays of one element
    else:
        delta = np.subtract(sample,old_mean)
        new_mean = old_mean + np.divide(delta,nsamples)
        new_M2 = old_M2 + delta*np.subtract(sample,new_mean)
        new_sample_variance = np.divide(new_M2,np.subtract(nsamples,1))
        if (compute_M3):
            new_M3 = old_M3 - 3.0*old_M2*np.divide(delta,nsamples) + np.divide(np.multiply((nsamples-1)*(nsamples-2),(delta**3)),(nsamples**2))
        else:
            new_M3 = old_M3 # we are not updating
        if (compute_M4):
            new_M4 = old_M4 - 4.0*old_M3*np.divide(delta,nsamples) + 6.0*old_M2*np.divide(delta,nsamples)**2 + np.multiply((nsamples-1)*(nsamples**2-3*nsamples+3),np.divide(delta**4,nsamples**3))
        else:
            new_M4 = old_M4 # we are not updating
    return new_mean,new_sample_variance,new_M2,new_M3,new_M4,nsamples


'''
auxiliary function of UpdateOnepassPowerSums of the StatisticalVariable class
'''
@ExaquteTask(returns=5)
def UpdateOnePassPowerSumsAux_Task(sample,old_S1,old_S2,old_S3,old_S3_absolute,old_S4,nsamples):
    if nsamples == 1:
        new_S1 = sample
        new_S2 = sample**2
        new_S3 = sample**3
        new_S3_absolute = np.abs(sample**3)
        new_S4 = sample**4
    else:
        new_S1 = old_S1 + sample
        new_S2 = old_S2 + sample**2
        new_S3 = old_S3 + sample**3
        new_S3_absolute = old_S3_absolute + np.abs(sample**3)
        new_S4 = old_S4 + sample**4
    return new_S1,new_S2,new_S3,new_S3_absolute,new_S4


'''
auxiliary function of UpdateHStatistics of the StatisticalVariable class
'''
@ExaquteTask(returns=4)
def ComputeHStatisticsAux_Task(S1_level,S2_level,S3_level,S4_level,number_samples_level):
    h1_level = 0
    h2_level = (number_samples_level*S2_level-S1_level**2) / ((number_samples_level-1)*number_samples_level)
    h3_level = (number_samples_level**2*S3_level-3*number_samples_level*S2_level*S1_level+2*S1_level**3) / \
        ((number_samples_level-2)*(number_samples_level-1)*number_samples_level)
    h4_level = ((-4*number_samples_level**2+8*number_samples_level-12)*S3_level*S1_level+ \
        (number_samples_level**3-2*number_samples_level**2+3*number_samples_level)*S4_level+ \
        6*number_samples_level*S2_level*S1_level**2+(9-6*number_samples_level)*S2_level**2-3*S1_level**4) / \
        ((number_samples_level-3)*(number_samples_level-2)*(number_samples_level-1)*number_samples_level)
    return h1_level,h2_level,h3_level,h4_level


'''
auxiliary function of ComputeSkewnessKurtosis of the StatisticalVariable class
'''
@ExaquteTask(returns=2)
def ComputeSkewnessKurtosisAux_Task(h2_level,h3_level,h4_level):
    skewness_level = h3_level / (np.sqrt(h2_level**3))
    kurtosis_level = h4_level / (h2_level**2)
    return skewness_level,kurtosis_level


'''
auxiliary function of ComputeSampleCentralMoments of the StatisticalVariable class
'''
@ExaquteTask(returns=5)
def ComputeSampleCentralMomentsAux_Task(sample,curr_mean,number_samples_level,compute_first_central_moment,compute_second_central_moment, \
    compute_third_central_moment,compute_third_absolute_central_moment,compute_fourth_central_moment, \
    first_central_moment,second_central_moment,third_central_moment,third_central_moment_absolute,fourth_central_moment):
    if (compute_first_central_moment):
        first_central_moment = first_central_moment + ((sample - curr_mean)**1) / number_samples_level
    if (compute_second_central_moment):
        second_central_moment = second_central_moment + ((sample - curr_mean)**2) / number_samples_level
    if (compute_third_central_moment):
        third_central_moment = third_central_moment + ((sample - curr_mean)**3) / number_samples_level
    if (compute_third_absolute_central_moment):
        third_central_moment_absolute = third_central_moment_absolute + (np.abs(sample - curr_mean)**3) / number_samples_level
    if (compute_fourth_central_moment):
        fourth_central_moment = fourth_central_moment + ((sample - curr_mean)**4) / number_samples_level
    return first_central_moment,second_central_moment,third_central_moment,third_central_moment_absolute,fourth_central_moment


class StatisticalVariable(object):
    '''The base class for statistical variables'''
    def __init__(self, number_levels):
        '''constructor of the class
        Keyword arguments:
        self : an instance of a class
        number_levels : number of levels
        '''

        '''values of the variable, organized per level'''
        self.values = []
        '''mean of the variable per each level'''
        self.mean = []
        '''sample variance of the variable per each level'''
        self.sample_variance = []
        '''moments of the variable per each level M_p   = n * \mu_p
                                                  \mu_p = p-th central moment
                                                  n     = number of values'''
        self.moment_2 = []
        self.moment_3 = []
        self.moment_4 = []
        '''set if M3 and M4 will be computed (M2 is mandatory to be computed because it is exploited in the mean evaluation)'''
        self.moment_3_to_compute = False
        self.moment_4_to_compute = False
        '''bias error of the variable'''
        self.bias_error = None
        '''statistical error of the variable'''
        self.statistical_error = None
        '''type of variable: scalar or field'''
        self.type = None
        '''number of samples of the variable'''
        self.number_samples = [0 for _ in range(number_levels+1)]
        '''power sums S_p = \sum_{i=1}^{n} Q(sample_i)**p, organized per level'''
        self.power_sum_1 = []
        self.power_sum_2 = []
        self.power_sum_3 = []
        self.power_sum_3_absolute = [] # S_p = \sum_{i=1}^{n} abs(Q(sample_i)**p)
        self.power_sum_4 = []
        '''sample central moments \mu_p = \sum_{i=1}^{n} (Q(sample_i)-mean_n)**p / n, organized per level'''
        self.sample_central_moment_1 = []
        self.sample_central_moment_2 = []
        self.sample_central_moment_3 = []
        self.sample_central_moment_3_absolute = [] # \mu_p = \sum_{i=1}^{n} abs((Q(sample_i)-mean_n)**p) / n
        self.sample_central_moment_4 = []
        self.sample_first_central_moment_to_compute = False
        self.sample_second_central_moment_to_compute = False
        self.sample_third_central_moment_to_compute = False
        self.sample_third_absolute_central_moment_to_compute = False
        self.sample_fourth_central_moment_to_compute = False
        '''h-statistics h_p, the unbiased central moment estimator with minimal variance, organized per level'''
        self.h_statistics_1 = []
        self.h_statistics_2 = []
        self.h_statistics_3 = []
        self.h_statistics_4 = []
        self.h_statistics_computed = False
        '''skewness of the variable per each level'''
        self.skewness = []
        '''kurtosis of the variable per each level'''
        self.kurtosis = []
        '''convergence criteria of the algorithm'''
        self.convergence_criteria = None

    '''
    function initializing variables of the Statistical Variable class given number of levels
    '''
    def InitializeStatisticalVariableVariables(self,number_levels):
        self.values = [[] for _ in range (number_levels)]
        self.mean = [[] for _ in range (number_levels)]
        self.moment_2 = [[] for _ in range (number_levels)]
        self.moment_3 = [[] for _ in range (number_levels)]
        self.moment_4 = [[] for _ in range (number_levels)]
        self.sample_variance = [[] for _ in range (number_levels)]
        self.power_sum_1 = [[] for _ in range (number_levels)]
        self.power_sum_2 = [[] for _ in range (number_levels)]
        self.power_sum_3 = [[] for _ in range (number_levels)]
        self.power_sum_3_absolute = [[] for _ in range (number_levels)]
        self.power_sum_4 = [[] for _ in range (number_levels)]
        self.h_statistics_1 = [[] for _ in range (number_levels)]
        self.h_statistics_2 = [[] for _ in range (number_levels)]
        self.h_statistics_3 = [[] for _ in range (number_levels)]
        self.h_statistics_4 = [[] for _ in range (number_levels)]
        self.skewness = [[] for _ in range (number_levels)]
        self.kurtosis = [[] for _ in range (number_levels)]
        self.sample_central_moment_1 = [[] for _ in range (number_levels)]
        self.sample_central_moment_2 = [[] for _ in range (number_levels)]
        self.sample_central_moment_3 = [[] for _ in range (number_levels)]
        self.sample_central_moment_3_absolute = [[] for _ in range (number_levels)]
        self.sample_central_moment_4 = [[] for _ in range (number_levels)]

    '''
    function updating moments and sample variance
    moments:
    M_p   = n * \mu_p
    \mu_p = p-th central moment
    n     = number of values
    sample variance:
    s_n^2 = M_{2,n} / (n-1)
    '''
    def UpdateOnePassMomentsVariance(self,level,i_sample):
        sample = self.values[level][i_sample]
        old_mean = self.mean[level]
        old_M2 = self.moment_2[level]
        old_M3 = self.moment_3[level]
        compute_M3 = self.moment_3_to_compute
        old_M4 = self.moment_4[level]
        compute_M4 = self.moment_4_to_compute
        number_samples_level = self.number_samples[level]
        new_mean,new_sample_variance,new_M2,new_M3,new_M4,number_samples_level = UpdateOnePassMomentsVarianceAux_Task(sample,old_mean,old_M2,old_M3,compute_M3,old_M4,compute_M4,number_samples_level)
        self.mean[level] = new_mean
        self.sample_variance[level] = new_sample_variance
        self.moment_2[level] = new_M2
        self.moment_3[level] = new_M3
        self.moment_4[level] = new_M4
        self.number_samples[level] = number_samples_level

    '''
    function updating the power sums S_p
    '''
    def UpdateOnePassPowerSums(self,level,i_sample):
        sample = self.values[level][i_sample]
        old_S1 = self.power_sum_1[level]
        old_S2 = self.power_sum_2[level]
        old_S3 = self.power_sum_3[level]
        old_S3_absolute = self.power_sum_3_absolute[level]
        old_S4 = self.power_sum_4[level]
        number_samples_level = self.number_samples[level]
        new_S1,new_S2,new_S3,new_S3_absolute,new_S4 = UpdateOnePassPowerSumsAux_Task(sample,old_S1,old_S2,old_S3,old_S3_absolute,old_S4,number_samples_level)
        self.power_sum_1[level] = new_S1
        self.power_sum_2[level] = new_S2
        self.power_sum_3[level] = new_S3
        self.power_sum_3_absolute[level] = new_S3_absolute
        self.power_sum_4[level] = new_S4

    '''
    function computing the h statistics h_p from the power sums
    '''
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

    '''
    function computing the central moments (and the absolute third central moment)
    '''
    def ComputeSampleCentralMoments(self,level,number_samples_level):
        '''local variables for mean and number of samples'''
        curr_mean = self.mean[level]
        '''initialize central moements'''
        first_central_moment = 0.0
        second_central_moment = 0.0
        third_central_moment = 0.0
        third_central_moment_absolute = 0.0
        fourth_central_moment = 0.0
        compute_first_central_moment = self.sample_first_central_moment_to_compute
        compute_second_central_moment = self.sample_second_central_moment_to_compute
        compute_third_central_moment = self.sample_third_central_moment_to_compute
        compute_third_absolute_central_moment = self.sample_third_absolute_central_moment_to_compute
        compute_fourth_central_moment = self.sample_fourth_central_moment_to_compute
        for i in range(0,number_samples_level):
            '''compute only the central moements we need, since it is expensive their computation at large number_samples_level'''
            sample = self.values[level][i]
            first_central_moment,second_central_moment,third_central_moment,third_central_moment_absolute,fourth_central_moment = \
                ComputeSampleCentralMomentsAux_Task(sample,curr_mean,number_samples_level,compute_first_central_moment, \
                compute_second_central_moment,compute_third_central_moment,compute_third_absolute_central_moment,compute_fourth_central_moment, \
                first_central_moment,second_central_moment,third_central_moment,third_central_moment_absolute,fourth_central_moment)
        self.sample_central_moment_1[level] = first_central_moment
        self.sample_central_moment_2[level] = second_central_moment
        self.sample_central_moment_3[level] = third_central_moment
        self.sample_central_moment_3_absolute[level] = third_central_moment_absolute
        self.sample_central_moment_4[level] = fourth_central_moment

    '''
    function computing the skewness and the kurtosis from the central moments (i.e. from the h statistics)
    skewness = \mu_3 / \sqrt(\mu_2^3)
    kurtosis = \mu_4 / \mu_2^2
    '''
    def ComputeSkewnessKurtosis(self,level):
        if (self.h_statistics_computed):
            h2_level = self.h_statistics_2[level]
            h3_level = self.h_statistics_3[level]
            h4_level = self.h_statistics_4[level]
            skewness_level,kurtosis_level =ComputeSkewnessKurtosisAux_Task(h2_level,h3_level,h4_level)
            self.skewness[level] = skewness_level
            self.kurtosis[level] = kurtosis_level
