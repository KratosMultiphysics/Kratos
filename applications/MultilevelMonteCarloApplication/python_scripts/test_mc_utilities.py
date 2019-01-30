from __future__ import absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import packages
import numpy as np
from math import *

# Import the StatisticalVariable class
from test_auxiliary_classes_utilities import StatisticalVariable

# Import exaqute
from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with pycompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
# from exaqute.ExaquteTaskLocal import *      # to execute with python3
'''
get_value_from_remote is the equivalent of compss_wait_on: a synchronization point
in future, when everything is integrated with the it4i team, importing exaqute.ExaquteTaskHyperLoom you can launch your code with their scheduler instead of BSC
'''

'''
This utility contains the functions to perform the Monte Carlo (CMLMC) algorithm

References:
M. Pisaroni, F. Nobile, P. Leyland; A Continuation Multi Level Monte Carlo (C-MLMC) method for uncertainty quantification in compressible inviscid aerodynamics; Computer Methods in Applied Mechanics and Engineering, vol 326, pp 20-50, 2017. DOI : 10.1016/j.cma.2017.07.030.
M. Pisaroni, S. Krumscheid, F. Nobile; Quantifying uncertain system outputs via the multilevel Monte Carlo method - Part I: Central moment estimation ;  available as MATHICSE technical report no. 23.2017
C. Bayer, H. Hoel, E. von Schwerin, R. Tempone; On NonAsymptotyc optimal stopping criteria in Monte Carlo simulations; avaiable at  SIAM Journal on Scientific Computing, 2014, Vol. 36, No. 2 : pp. A869-A885
'''


'''
auxiliary function of AddResults of the MonteCarlo class
'''
def AddResultsAux_Task(simulation_results,level):
    if (level == 0):
        QoI_value = simulation_results
    else:
        raise Exception("level not equal to 0, in MC we should have only level zero")
    return QoI_value


'''
auxiliary function of CheckConvergence of the MonteCarlo class
'''
def CheckConvergenceAux_Task(curr_number_samples,curr_mean,curr_sample_variance,curr_h2,curr_h3,curr_sample_central_moment_3_absolute,curr_h4,curr_tol,curr_delta,convergence_criteria):
    convergence_boolean = False
    if(convergence_criteria == "MC_sample_variance_sequential_stopping_rule"):
        '''define local variables'''
        curr_coefficient_to_compute_convergence = np.sqrt(curr_number_samples) * curr_tol / np.sqrt(curr_sample_variance)
        '''evaluate probability of failure'''
        main_contribute = 2*(1-_ComputeCDFStandardNormalDistribution(curr_coefficient_to_compute_convergence))
        if(main_contribute < curr_delta):
            convergence_boolean = True
    elif(convergence_criteria == "MC_higher_moments_sequential_stopping_rule"):
        '''define local variables'''
        curr_second_sample_moment = np.sqrt(curr_h2)
        curr_third_sample_moment_absolute = curr_sample_central_moment_3_absolute / (curr_second_sample_moment**3)
        curr_third_sample_moment = curr_h3 / (curr_second_sample_moment**3)
        curr_fourth_sample_moment = (curr_h4 / (curr_second_sample_moment**4)) - 3
        curr_coefficient_to_compute_convergence = np.sqrt(curr_number_samples) * curr_tol / np.sqrt(curr_second_sample_moment)
        '''evaluate probability of failure and penalty term'''
        main_contribute = 2 * (1 - _ComputeCDFStandardNormalDistribution(curr_coefficient_to_compute_convergence))
        penalty_contribute = 2 * np.minimum(4 * (2/(curr_number_samples - 1) + (curr_fourth_sample_moment / curr_number_samples)), 1) * \
            _ComputeBoundFunction(curr_coefficient_to_compute_convergence,curr_third_sample_moment_absolute) / np.sqrt(curr_number_samples) + \
            (1 - np.minimum(4 * (2/(curr_number_samples - 1) + (curr_fourth_sample_moment / curr_number_samples)), 1)) * \
            np.abs((curr_number_samples * curr_tol**2 / (curr_second_sample_moment**2)) - 1) * \
            np.exp(- curr_number_samples * (curr_tol**2) / (curr_second_sample_moment**2)) * np.abs(curr_third_sample_moment) / \
            (3 * np.sqrt(2 * np.pi * curr_number_samples))
        if (main_contribute + penalty_contribute < curr_delta):
            convergence_boolean = True
    else:
        convergence_boolean = False
    return convergence_boolean


class MonteCarlo(object):
    '''The base class for the MonteCarlo-classes'''
    def __init__(self,custom_settings):
        '''constructor of the MonteCarlo-Object
        Keyword arguments:
        self     : an instance of a class
        settings : the settings of the Monte Carlo simulation
        '''

        '''
        Kratos Parameters object containing the default setings of the Continuation MC simulation
        tolerance : Tolerance final
        cphi : Confidence on tolerance
        batch_size : Number of samples per batch size
        '''
        default_settings = KratosMultiphysics.Parameters("""
        {
            "tolerance" : 1e-1,
            "cphi" : 1e-1,
            "batch_size" : 50,
            "convergence_criteria" : "MC_higher_moments_sequential_stopping_rule"
        }
        """)
        self.settings = custom_settings
        '''validate and assign default parameters'''
        self.settings.ValidateAndAssignDefaults(default_settings)
        '''convergence : boolean variable defining if MLMC algorithm is converged'''
        self.convergence = False
        '''current_number_levels : number of levels of MC by default == 0 (we only have level 0)'''
        self.current_number_levels = 0
        '''theta_i : splitting parameter \in (0,1), this affects bias and statistical error in the computation of the total error'''
        self.theta_i = None
        '''mean_mlmc_QoI : MLMC estimator for the mean value of the Quantity of Interest'''
        self.mean_mc_QoI = None
        '''TErr : total error of MC algorithm, the sum of bias and statistical error is an overestmation of the real total error
                  TErr := \abs(E^MC[QoI] - E[QoI])'''
        self.TErr = None
        '''QoI : Quantity of Interest of the considered problem'''
        self.QoI = StatisticalVariable(self.current_number_levels)
        '''initialize properly all the variables of the StatisticalVariable class: MC has only one level (level 0)'''
        self.QoI.values = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.mean = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.second_moment = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.sample_variance = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.power_sum_1 = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.power_sum_2 = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.power_sum_3 = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.power_sum_3_absolute = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.power_sum_4 = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.h_statistics_1 = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.h_statistics_2 = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.h_statistics_3 = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.h_statistics_4 = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.skewness = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.kurtosis = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.sample_central_moment_1 = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.sample_central_moment_2 = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.sample_central_moment_3 = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.sample_central_moment_3_absolute = [[] for _ in range (self.current_number_levels+1)]
        self.QoI.sample_central_moment_4 = [[] for _ in range (self.current_number_levels+1)]
        '''number_samples : total number of samples of current iteration'''
        self.number_samples = [self.settings["batch_size"].GetInt() for _ in range (self.current_number_levels+1)]
        '''difference_number_samples : difference between number of samples of current and previous iterations'''
        self.difference_number_samples = [self.settings["batch_size"].GetInt() for _ in range (self.current_number_levels+1)]
        '''previous_number_samples : total number of samples of previous iteration'''
        self.previous_number_samples = [0 for _ in range (self.current_number_levels+1)]
        '''iteration counter'''
        self.iteration_counter = 0
        '''set convergence criteria'''
        self.SetConvergenceCriteria(self.settings["convergence_criteria"].GetString())
        if (self.convergence_criteria != "MC_sample_variance_sequential_stopping_rule" and self.convergence_criteria != "MC_higher_moments_sequential_stopping_rule"):
            raise Exception ("The selected convergence criteria is not yet implemented, plese select one of the following: \n i)  MC_sample_variance_sequential_stopping_rule \n ii) MC_higher_moments_sequential_stopping_rule")

    '''
    function checking if the MC algorithm is conveged, with respect to the selected convergence criteria
    '''
    def CheckConvergence(self,level):
        curr_number_samples = self.QoI.number_samples[level]
        curr_mean = self.QoI.mean[level]
        curr_sample_variance = self.QoI.sample_variance[level]
        curr_h2 = self.QoI.h_statistics_2[level]
        curr_h3 = self.QoI.h_statistics_3[level]
        curr_sample_central_moment_3_absolute = self.QoI.sample_central_moment_3_absolute[level]
        curr_h4 = self.QoI.h_statistics_4[level]
        curr_tol = self.settings["tolerance"].GetDouble()
        curr_delta = self.settings["cphi"].GetDouble()
        convergence_criteria = self.convergence_criteria

        convergence_boolean = CheckConvergenceAux_Task(curr_number_samples,curr_mean,curr_sample_variance,curr_h2,\
            curr_h3,curr_sample_central_moment_3_absolute,curr_h4,curr_tol,curr_delta,convergence_criteria)
        self.convergence = convergence_boolean

    '''
    function adding QoI values to the corresponding level
    '''
    def AddResults(self,simulation_results):
        '''simulation_results = [MultilevelMonteCarloResults class, level (integer type, not compss.future)]'''
        level = 0
        QoI_value = AddResultsAux_Task(simulation_results,level)
        '''update values of QoI'''
        self.QoI.values[level] = np.append(self.QoI.values[level],QoI_value)

    '''
    function initializing the MC phase
    usage: It is designed to be called BEFORE the MC solution step
    '''
    def InitializeMCPhase(self):
        level = 0
        '''update iteration counter'''
        self.iteration_counter = self.iteration_counter + 1
        '''update number of samples'''
        if (self.iteration_counter > 1):
            self.previous_number_samples[level] = self.number_samples[level]
            self.number_samples[level] = self.number_samples[level] * 2 + self.previous_number_samples[level]
            self.difference_number_samples[level] = self.number_samples[level] - self.previous_number_samples[level]
        else:
            pass

    '''
    function finalizing the MC phase
    usage: It is designed to be called AFTER the MC solution step
    '''
    def FinalizeMCPhase(self):
        level = 0 # level = 0 (Monte Carlo algorithm)
        '''update statistics of the QoI'''
        for i_sample in range(self.previous_number_samples[level],self.number_samples[level]):
            self.QoI.UpdateOnepassMeanVariance(level,i_sample)
            self.QoI.UpdateOnePassPowerSums(level,i_sample)
        '''compute the central moments we can't derive from the unbiased h statistics
        we compute expensively the absolute central moment because we can't retrieve it from the h statistics'''
        self.QoI.sample_third_absolute_central_moment_to_compute = True
        self.QoI.ComputeSampleCentralMoments(level,self.number_samples[level]) # give as input the number of samples computed up to this point
        self.QoI.ComputeHStatistics(level)
        self.QoI.ComputeSkewnessKurtosis(level)
        self.CheckConvergence(level)
        '''synchronization point needed to launch new tasks, if needed
        put as in the end as possible the synchronization point'''
        self.convergence = get_value_from_remote(self.convergence)

    '''
    function printing informations about initializing MLMC phase
    '''
    def ScreeningInfoInitializeMCPhase(self):
        print("\n","#"*50," MC iter =  ",self.iteration_counter,"#"*50,"\n")

    '''
    function printing informations about finalizing MC phase
    '''
    def ScreeningInfoFinalizeMCPhase(self):
        # print("values computed of QoI = ",self.QoI.values)
        print("samples computed in this iteration",self.difference_number_samples)
        print("current number of samples = ",self.number_samples)
        print("previous number of samples = ",self.previous_number_samples)
        print("monte carlo mean and variance QoI estimators = ",self.QoI.mean,self.QoI.sample_variance)
        print("convergence = ",self.convergence)

    '''
    function setting the convergence criteria the algorithm will exploit
    '''
    def SetConvergenceCriteria(self,convergence_string_name):
        self.convergence_criteria = convergence_string_name


'''
auxiliary function of CheckConvergence for the MC_higher_moments_sequential_stopping_rule criteria
'''
def _ComputeBoundFunction(x,beta):
    return np.minimum(0.3328 * (beta + 0.429), 18.1139 * beta / (1 + (np.abs(x)**3)))


'''
function computing the cumulative distribution function for the standard normal distribution
'''
def _ComputeCDFStandardNormalDistribution(x):
    'cumulative distribution function (CDF) for the standard normal distribution'
    return (1.0 + erf(x / sqrt(2.0))) / 2.0