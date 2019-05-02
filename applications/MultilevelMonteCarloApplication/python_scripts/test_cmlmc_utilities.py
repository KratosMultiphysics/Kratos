from __future__ import absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import packages
import numpy as np
import copy

# Import the StatisticalVariable class
from test_auxiliary_classes_utilities import StatisticalVariable


'''
This utility contains the functions to perform the Continuation Multilevel Monte Carlo (CMLMC) algorithm

References:
M. Pisaroni, F. Nobile, P. Leyland; A Continuation Multi Level Monte Carlo (C-MLMC) method for uncertainty quantification in compressible inviscid aerodynamics; Computer Methods in Applied Mechanics and Engineering, vol 326, pp 20-50, 2017. DOI : 10.1016/j.cma.2017.07.030.
M. Pisaroni, S. Krumscheid, F. Nobile; Quantifying uncertain system outputs via the multilevel Monte Carlo method - Part I: Central moment estimation ;  available as MATHICSE technical report no. 23.2017
'''


'''
TODO: create a class for the values ValueClass in order to directly modify the values inside functions
e.g. mean = [1,2,4,5]  ---> mean = [Value, Value, Value, Value]
     to modify this inside a task I need it to be an object
TODO: insert distinction between scalar and field Quantity of Interests
'''


'''
auxiliary function of AddResults of the MultilevelMonteCarlo class
'''
def AddResultsAux_Task(simulation_results,level):
    if (level == 0):
        difference_QoI_value = simulation_results.QoI[level]
    else:
        difference_QoI_value = simulation_results.QoI[level] - simulation_results.QoI[level - 1]
    return difference_QoI_value,simulation_results.time_ML[level]


'''
auxiliary function finalizing the screening phase and the MLMC phase of the MultilevelMonteCarlo class
'''
def FinalizePhase_Task(ConstructorCallback,aux_settings_serialized,aux_mesh_parameters,\
aux_current_number_levels,aux_current_iteration,aux_number_samples,*args):
    '''retrieve lists'''
    args = list(args)
    first_occ = args.index("%%%")
    difference_QoI_mean = args[:first_occ]
    second_occ = args[first_occ + 1:].index("%%%")
    difference_QoI_sample_variance = args[first_occ + 1:first_occ + 1 + second_occ]
    time_ML_mean = args[first_occ + 1 + second_occ + 1:]
    '''load settings (Kratos Parameters object) from SteramSerializer Kratos object'''
    aux_settings = KratosMultiphysics.Parameters()
    aux_settings_serialized.Load("ParametersSerialization",aux_settings)
    '''create an auxiliary object auxiliary_MLMC_object equivalent to the current one of the problem
    construct the class with the aux_settings settings passed from outside'''
    auxiliary_MLMC_object = ConstructorCallback(aux_settings)
    auxiliary_MLMC_object.difference_QoI.mean = difference_QoI_mean
    auxiliary_MLMC_object.difference_QoI.sample_variance = difference_QoI_sample_variance
    auxiliary_MLMC_object.time_ML.mean = time_ML_mean
    auxiliary_MLMC_object.mesh_parameters = aux_mesh_parameters
    auxiliary_MLMC_object.current_number_levels = aux_current_number_levels
    auxiliary_MLMC_object.current_iteration = aux_current_iteration
    auxiliary_MLMC_object.number_samples = aux_number_samples
    '''compute the functions needed to finalize the screening phase or the MLMC phase'''
    if (auxiliary_MLMC_object.current_iteration == 0): # screening phase
        '''compute parameters by least square fit'''
        auxiliary_MLMC_object.ComputeRatesLS()
        '''compute Bayesian variance V^c[Y_l]'''
        auxiliary_MLMC_object.EstimateBayesianVariance(auxiliary_MLMC_object.current_number_levels)
    else: # MLMC phase
        '''compute estimator MLMC mean QoI'''
        auxiliary_MLMC_object.ComputeMeanMLMCQoI()
        '''compute parameters by least square fit'''
        auxiliary_MLMC_object.ComputeRatesLS()
        '''compute Bayesian variance V^c[Y_l]'''
        auxiliary_MLMC_object.EstimateBayesianVariance(auxiliary_MLMC_object.current_number_levels)
        '''compute total error of MLMC simulation'''
        auxiliary_MLMC_object.ComputeTotalErrorMLMC()
    return auxiliary_MLMC_object.rates_error,\
    auxiliary_MLMC_object.BayesianVariance,auxiliary_MLMC_object.mean_mlmc_QoI,\
    auxiliary_MLMC_object.TErr,auxiliary_MLMC_object.number_samples


class MultilevelMonteCarlo(object):
    '''The base class for the MultilevelMonteCarlo-classes'''
    def __init__(self,custom_settings):
        '''constructor of the MultilevelMonteCarlo-Object
        Keyword arguments:
        self     : an instance of a class
        settings : the settings of the Multilevel Monte Carlo simulation
        '''

        '''
        Kratos Parameters object containing the default setings of the Continuation MLMC simulation
        tol0 : Tolerance iter 0
        tolF : Tolerance final
        cphi : Confidence on tolerance
        number_samples_screening : Number of samples for iter 0
        L0   : Number of levels for iter 0
        Lmax : Maximum number of levels
        mesh_refinement_coefficient : coefficient of mesh refinement
        initial_mesh_size : size of first mesh considered
        minimum_add_level : minimum number of samples to add if at least one is going to be added
        k0   : Certainty Parameter 0 rates (confidence in the variance models)
        k1   : Certainty Parameter 1 rates (confidence in the weak convergence model)
        r1   : Cost increase first iterations C-MLMC
        r2   : Cost increase final iterations C-MLMC
        '''
        default_settings = KratosMultiphysics.Parameters("""
        {
            "k0" : 0.1,
            "k1" : 0.1,
            "r1" : 1.25,
            "r2" : 1.15,
            "tol0" : 0.25,
            "tolF" : 0.1,
            "cphi" : 1.0,
            "number_samples_screening" : 25,
            "Lscreening" : 2,
            "Lmax" : 4,
            "mesh_refinement_coefficient" : 2,
            "initial_mesh_size" : 0.5,
            "minimum_add_level" : 6.0,
            "splitting_parameter_max" : 0.9,
            "splitting_parameter_min" : 0.1
        }
        """)
        self.settings = custom_settings
        '''warning if initial_mesh_size parameter not set by the user'''
        if not (self.settings.Has("initial_mesh_size")):
            print("\n ######## WARNING : initial_mesh_size parameter not set ---> using defalut value 0.5 ########\n")
        '''validate and assign default parameters'''
        self.settings.ValidateAndAssignDefaults(default_settings)
        '''current_number_levels : number of levels of current iteration'''
        self.current_number_levels = self.settings["Lscreening"].GetInt()
        '''previous_number_levels : number of levels of previous iteration'''
        self.previous_number_levels = None
        '''number_samples : total number of samples of current iteration'''
        self.number_samples = [self.settings["number_samples_screening"].GetInt() for _ in range (self.settings["Lscreening"].GetInt()+1)]
        '''difference_number_samples : difference between number of samples of current and previous iterations'''
        self.difference_number_samples = None
        '''previous_number_samples : total number of samples of previous iteration'''
        self.previous_number_samples = None
        '''rates_error : dictionary containing the values of the parameters
        calpha : coefficient of the function maximizing bias
        alpha  : exponent of the function maximizing bias
        cbeta  : coefficient of the function maximizing statistical error
        beta   : exponent of the function maximizing statistical error
        cgamma : coefficient of the function maximizing cost
        gamma  : exponent of the function maximizing cost
        '''
        self.rates_error = {"calpha":None, "alpha":None, "cbeta":None, "beta":None, "cgamma":None, "gamma":None}
        '''mesh_parameters : reciprocal of minimal mesh size'''
        self.mesh_parameters = []
        '''size_mesh : minimal mesh size'''
        self.sizes_mesh = []
        '''BayesianVariance : Bayesian variance'''
        self.BayesianVariance = []
        '''number_iterations : theoretical number of iterations the MLMC algorithm will perform'''
        self.number_iterations_iE = None
        '''convergence : boolean variable defining if MLMC algorithm is converged'''
        self.convergence = False
        '''current_iteration : current iteration of MLMC algorithm'''
        self.current_iteration = 0
        '''tolerance_i : tolerance of i^th-iteration of MLMC algorithm'''
        self.tolerance_i = None
        '''theta_i : splitting parameter \in (0,1), this affects bias and statistical error in the computation of the total error'''
        self.theta_i = None
        '''mean_mlmc_QoI : MLMC estimator for the mean value of the Quantity of Interest'''
        self.mean_mlmc_QoI = None
        '''TErr : total error of MLMC algorithm, the sum of bias and statistical error is an overestmation of the real total error
                  TErr := \abs(E^MLMC[QoI] - E[QoI])'''
        self.TErr = None
        '''compute mesh parameter for each mesh'''
        self.ComputeMeshParameters()

        '''difference_QoI : Quantity of Interest of the considered problem organized per levels
                            difference_QoI.values := Y_l = QoI_M_l - Q_M_l-1'''
        self.difference_QoI = StatisticalVariable(self.current_number_levels)
        self.difference_QoI.values = [[] for _ in range (self.settings["Lscreening"].GetInt()+1)] # list containing Y_{l}^{i} = Q_{m_l} - Q_{m_{l-1}}
        self.difference_QoI.type = "scalar"
        '''time_ML : time to perform a single MLMC simulation (i.e. one value of difference_QoI.values) organized per levels'''
        self.time_ML = StatisticalVariable(self.current_number_levels)
        self.time_ML.values = [[] for _ in range (self.settings["Lscreening"].GetInt()+1)] # list containing the time to compute the level=l simulations

        '''########################################################################
        # observation: levels start from level 0                                  #
        #              length arrays and lists starts from 1                      #
        # thus there is a difference of 1 between length lists and levels         #
        # e.g. self.current_level = len(self.number_samples) - 1                  #
        #      or                                                                 #
        #      self.current_level = len(self.difference_QoI.values) - 1           #
        ########################################################################'''

    '''
    function finalizing the screening phase of MLMC algorithm
    Usage: It is designed to be called ONCE, AFTER the screening phase
    '''
    def FinalizeScreeningPhase(self):
        '''prepare lists'''
        self.difference_QoI.mean = [[] for _ in range (self.settings["Lscreening"].GetInt()+1)]
        self.difference_QoI.sample_variance = [[] for _ in range (self.settings["Lscreening"].GetInt()+1)]
        self.difference_QoI.second_moment = [[] for _ in range (self.settings["Lscreening"].GetInt()+1)]
        self.time_ML.mean = [[] for _ in range (self.settings["Lscreening"].GetInt()+1)]
        self.time_ML.sample_variance = [[] for _ in range (self.settings["Lscreening"].GetInt()+1)]
        self.time_ML.second_moment = [[] for _ in range (self.settings["Lscreening"].GetInt()+1)]
        '''compute mean, sample variance and second moment for difference QoI and time ML'''
        for level in range (self.current_number_levels+1):
            for i_sample in range(self.number_samples[level]):
                self.difference_QoI.UpdateOnepassMeanVariance(level,i_sample)
                self.time_ML.UpdateOnepassMeanVariance(level,i_sample)
        '''compute i_E, number of iterations of Multilevel Monte Carlo algorithm'''
        self.ComputeNumberIterationsMLMC()
        '''the following lines represent the functions we execute in FinalizePhase_Task
        self.ComputeRatesLS()
        self.EstimateBayesianVariance(self.current_number_levels)'''
        '''store lists in synchro_element to execute FinalizePhase_Task'''
        synchro_elements = [x for x in self.difference_QoI.mean]
        synchro_elements.extend(["%%%"])
        synchro_elements.extend(self.difference_QoI.sample_variance)
        synchro_elements.extend(["%%%"])
        synchro_elements.extend(self.time_ML.mean)
        '''create a StreamSerializer Kratos object containing the problem settings'''
        serial_settings = KratosMultiphysics.StreamSerializer()
        serial_settings.Save("ParametersSerialization", self.settings)
        '''compute remaining MLMC finalize process operations
        observation: we are passing self.settings and we will exploit it to construct the class'''
        self.rates_error,self.BayesianVariance,self.mean_mlmc_QoI,\
        self.TErr,self.number_samples\
        = FinalizePhase_Task(MultilevelMonteCarlo,
        serial_settings,self.mesh_parameters,self.current_number_levels,\
        self.current_iteration,self.number_samples,*synchro_elements)
        '''start first iteration, we enter in the MLMC algorithm'''
        self.current_iteration = 1

    '''
    function performing all the required operations that should be executed
    (for each MLMC iteration) BEFORE the MLMC solution step
    '''
    def InitializeMLMCPhase(self):
        '''compute tolerance for the i^th iteration'''
        self.ComputeTolerancei()
        '''Compute Optimal Number of Levels for iteration i L_i'''
        self.ComputeLevels()
        '''compute theta splitting parameter according to the current_number_levels and tolerance_i'''
        self.ComputeTheta(self.current_number_levels)
        '''compute number of samples according to BayesianVariance and theta_i parameters'''
        self.ComputeNumberSamples()
        '''prepare lists'''
        for _ in range (self.current_number_levels - self.previous_number_levels): # append a list for the new level
            self.difference_QoI.values.append([])
            self.time_ML.values.append([])

    '''function performing all the required operations that should be executed
    (for each MLMC iteration) AFTER the MLMC solution step'''
    def FinalizeMLMCPhase(self):
        '''prepare lists'''
        for _ in range (self.current_number_levels - self.previous_number_levels): # append a list for the new level
            self.difference_QoI.mean.append([])
            self.difference_QoI.sample_variance.append([])
            self.difference_QoI.second_moment.append([])
            self.difference_QoI.number_samples.append(0)
            self.time_ML.mean.append([])
            self.time_ML.sample_variance.append([])
            self.time_ML.second_moment.append([])
            self.time_ML.number_samples.append(0)
        '''compute mean, sample variance and second moment for difference QoI and time ML'''
        for level in range (self.current_number_levels+1):
            for i_sample in range(self.previous_number_samples[level],self.number_samples[level]):
                self.difference_QoI.UpdateOnepassMeanVariance(level,i_sample)
                self.time_ML.UpdateOnepassMeanVariance(level,i_sample)
        '''update number of levels'''
        self.previous_number_levels = self.current_number_levels
        '''the following commented lines represent the functions we execute in FinalizePhase_Task
        self.ComputeMeanMLMCQoI()
        self.ComputeRatesLS()
        self.EstimateBayesianVariance(self.current_number_levels)
        self.ComputeTotalErrorMLMC()'''
        '''store lists in synchro_element to execute FinalizePhase_Task'''
        synchro_elements = [x for x in self.difference_QoI.mean]
        synchro_elements.extend(["%%%"])
        synchro_elements.extend(self.difference_QoI.sample_variance)
        synchro_elements.extend(["%%%"])
        synchro_elements.extend(self.time_ML.mean)
        '''create a StreamSerializer Kratos object containing the problem settings'''
        serial_settings = KratosMultiphysics.StreamSerializer()
        serial_settings.Save("ParametersSerialization", self.settings)
        '''compute remaining MLMC finalize process operations
        observation: we are passing self.settings and we will exploit it to construct the class'''
        self.rates_error,self.BayesianVariance,self.mean_mlmc_QoI,\
        self.TErr,self.number_samples\
        = FinalizePhase_Task(MultilevelMonteCarlo,
        serial_settings,self.mesh_parameters,self.current_number_levels,\
        self.current_iteration,self.number_samples,*synchro_elements)

        '''check convergence
        convergence reached if: i) current_iteration >= number_iterations_iE
                               ii) TErr < tolerance_i
        if not update current_iteration'''
        if (self.current_iteration >= self.number_iterations_iE) and (self.TErr < self.tolerance_i):
            self.convergence = True
        else:
            self.current_iteration = self.current_iteration + 1

    '''
    function printing informations about screening phase
    '''
    def ScreeningInfoScreeningPhase(self):
        print("\n","#"*50," SCREENING PHASE ","#"*50,"\n")
        # print("values computed of QoI = ",self.difference_QoI.values)
        # print("values computed time_ML",self.time_ML.values)
        print("mean and variance difference_QoI = ",self.difference_QoI.mean,self.difference_QoI.sample_variance)
        print("mean time_ML",self.time_ML.mean)
        print("rates coefficient = ",self.rates_error)
        print("estimated Bayesian variance = ",self.BayesianVariance)
        print("minimum number of MLMC iterations = ",self.number_iterations_iE)

    '''
    function printing informations about initializing MLMC phase
    '''
    def ScreeningInfoInitializeMLMCPhase(self):
        print("\n","#"*50," MLMC iter =  ",self.current_iteration,"#"*50,"\n")
        print("current tolerance = ",self.tolerance_i)
        print("updated estimated Bayesian Variance initialize phase = ",self.BayesianVariance)
        print("current number of levels = ",self.current_number_levels)
        print("previous number of levels = ",self.previous_number_levels)
        print("current splitting parameter = ",self.theta_i)
        print("difference number of samples = ",self.difference_number_samples)
        print("previous number of samples = ",self.previous_number_samples)

    '''
    function printing informations about finalizing MLMC phase
    '''
    def ScreeningInfoFinalizeMLMCPhase(self):
        # print("values computed of QoI = ",self.difference_QoI.values)
        # print("values computed time_ML",self.time_ML.values)
        print("current number of samples",self.number_samples)
        print("mean and variance difference_QoI = ",self.difference_QoI.mean,self.difference_QoI.sample_variance)
        print("mean time_ML",self.time_ML.mean)
        print("rates coefficient = ",self.rates_error)
        print("estimated Bayesian variance = ",self.BayesianVariance)
        print("multilevel monte carlo mean estimator = ",self.mean_mlmc_QoI)
        print("TErr = bias + statistical error = ",self.TErr)

    '''
    function adding QoI and MLMC time values to the corresponding level
    '''
    def AddResults(self,simulation_results):
        '''simulation_results = [MultilevelMonteCarloResults class, level (integer type, not compss.future)]'''
        simulation_results_class = simulation_results[0]
        level = simulation_results[1]
        for lev in range (level+1):
            difference_QoI_value, time_ML_value = AddResultsAux_Task(simulation_results_class,lev)
            '''update values of difference QoI and time ML per each level'''
            self.difference_QoI.values[lev] = np.append(self.difference_QoI.values[lev], difference_QoI_value)
            self.time_ML.values[lev] = np.append(self.time_ML.values[lev],time_ML_value)
            '''update number of samples'''
            if (lev != level):
                self.number_samples[lev] = self.number_samples[lev] + 1

    '''
    function giving as output the mesh discretization parameter
    the mesh parameter is the reciprocal of the minimum mesh size of the grid
    h_lev=h_0*M^(-lev)
    '''
    def ComputeMeshParameters(self):
        h0 = self.settings["initial_mesh_size"].GetDouble()
        M  = self.settings["mesh_refinement_coefficient"].GetInt()
        for level in range(self.settings["Lmax"].GetInt()+1):
            h_current_level = h0 * M**(-level)
            mesh_parameter_current_level = h_current_level**(-1)
            self.sizes_mesh.append(h_current_level)
            self.mesh_parameters.append(mesh_parameter_current_level)

    '''
    function computing the problem parameters P=[calpha,alpha,cbeta,beta,cgamma,gamma] using least squares fit
    we consider level > 0 to compute calpha,alpha,cbeta,beta for robustness reasons [see PNL17 for details]
    '''
    def ComputeRatesLS(self):
        bias_ratesLS = np.abs(self.difference_QoI.mean)
        variance_ratesLS = self.difference_QoI.sample_variance
        cost_ML_ratesLS = self.time_ML.mean
        mesh_param_ratesLS = self.mesh_parameters[0:self.current_number_levels+1]
        '''mean - alpha
        linear fit'''
        pa = np.polyfit(np.log2(mesh_param_ratesLS[1::]),np.log2(bias_ratesLS[1::]),1)
        alpha   = -pa[0]
        C1      = 2**pa[1]
        '''variance - beta
        linear fit'''
        pb          = np.polyfit(np.log2(mesh_param_ratesLS[1::]),np.log2(variance_ratesLS[1::]),1)
        beta        = -pb[0]
        C2          = 2**pb[1]
        '''cost of computation - gamma
        linear fit'''
        pg          = np.polyfit(np.log2(mesh_param_ratesLS),np.log2(cost_ML_ratesLS),1)
        gamma       = pg[0]
        C3          = 2**pg[1]
        '''update the rates error dictionary'''
        self.rates_error["calpha"] = C1
        self.rates_error["alpha"] = alpha
        self.rates_error["cbeta"] = C2
        self.rates_error["beta"] = beta
        self.rates_error["cgamma"] = C3
        self.rates_error["gamma"] = gamma
        del(bias_ratesLS,variance_ratesLS,cost_ML_ratesLS,mesh_param_ratesLS,C1,C2,C3,alpha,beta,gamma)

    '''
    function performing the Bayesian update of the variance
    using samples generated on all levels in order to locally improve the estimation of variance[difference_QoI]
    '''
    def EstimateBayesianVariance(self,levels): # need to keep levels because in ComputeLevels I use the maximum number of levels
        '''use local variables'''
        k0 = self.settings["k0"].GetDouble()
        k1 = self.settings["k1"].GetDouble()
        calfa = self.rates_error["calpha"]
        alfa  = self.rates_error["alpha"]
        cbeta = self.rates_error["cbeta"]
        beta  = self.rates_error["beta"]
        mesh_param = self.mesh_parameters
        '''use local variables, in order to not modify the global variables'''
        mean_local = copy.copy(self.difference_QoI.mean)
        variance_local = copy.copy(self.difference_QoI.sample_variance)
        nsam_local = copy.copy(self.number_samples)
        '''append null values to evaluate the Bayesian variance for all levels'''
        if len(mean_local) < (levels+1):
            for _ in range (0,(levels+1)-len(mean_local)):
                mean_local.append(0.0)
        if len(variance_local) < (levels+1):
            for _ in range (0,(levels+1)-len(variance_local)):
                variance_local.append(0.0)
        if len(nsam_local) < (levels+1):
            for _ in range ((levels+1)-len(nsam_local)):
                nsam_local.append(0)
        '''estimate the Bayesian variance'''
        BayesianVariance = []
        for level in range (0, (levels+1)):
            mu = calfa*mesh_param[level]**(-alfa)
            lam = (1/cbeta)*mesh_param[level]**(beta)
            G1_l = 0.5 + np.multiply(k1,lam) + np.divide(nsam_local[level],2.0)
            G2_l = k1 + (nsam_local[level]-1)*0.5*variance_local[level] + k0*nsam_local[level]*((mean_local[level]-mu)**2)/(2.0*(k0+nsam_local[level]))
            BayesianVariance.append(np.divide(G2_l,G1_l-0.5))
        self.BayesianVariance = BayesianVariance
        del(BayesianVariance)

    '''
    function computing the minimum number of iterations of MLMC algorithm
    iteration < number_iterations_iE : iterations needed to obtain increasingly accurate estimates of the
                                       problem dependent parameters P=[calpha,alpha,cbeta,beta,cgamma,gamma]
    iteration > number_iterations_iE : iterations preventing redundant computations due to fluctuations
                                       in the estimate of P=[calpha,alpha,cbeta,beta,cgamma,gamma]
                                       by solving the problem for a slightly smaller tolerance than the desired one
    '''
    def ComputeNumberIterationsMLMC(self):
        tolF = self.settings["tolF"].GetDouble()
        tol0 = self.settings["tol0"].GetDouble()
        r2 = self.settings["r2"].GetDouble()
        r1 = self.settings["r1"].GetDouble()
        self.number_iterations_iE = np.floor((-np.log(tolF)+np.log(r2)+np.log(tol0))/(np.log(r1)))

    '''
    function computing the tolerance for the i^th iteration
    '''
    def ComputeTolerancei(self):
        tolF = self.settings["tolF"].GetDouble()
        r2 = self.settings["r2"].GetDouble()
        r1 = self.settings["r1"].GetDouble()
        iE = self.number_iterations_iE
        iter_def = self.current_iteration
        if iter_def < iE:
            tol = (r1**(iE-iter_def) * r2**(-1))*tolF
        elif iter_def > iE:
            tol = (r2**(iE-iter_def) * r2**(-1))*tolF
        else:
            tol = tolF
        self.tolerance_i = tol

    '''
    function computing the number of levels for i^th iteration of the algorithm
    '''
    def ComputeLevels(self):
        tol = self.tolerance_i
        Wmin   = 1e10 # high cost to compare with all the level costs (needs to be a high value)
        current_number_levels = self.current_number_levels
        cgamma = self.rates_error["cgamma"]
        gamma  = self.rates_error["gamma"]
        calpha = self.rates_error["calpha"]
        alpha = self.rates_error["alpha"]
        cphi = self.settings["cphi"].GetDouble()
        mesh_parameters_all_levels = self.mesh_parameters
        Lmax = self.settings["Lmax"].GetInt()
        Lmin = self.current_number_levels
        '''prepare mesh_parameters_all_levels and BayesianVariance to have both length = Lmax + 1'''
        if len(self.BayesianVariance) < (Lmax+1):
            self.EstimateBayesianVariance(Lmax)
        BayesianVariance = self.BayesianVariance
        model_cost = np.multiply(cgamma,np.power(mesh_parameters_all_levels,gamma)) # model_cost has length = Lmax + 1
        '''observation: not mandatory to increase the number of levels, we may continue using the number of levels of the previous MLMC iteration'''
        for lev in range(Lmin, Lmax+1):
            # theta_i = 1.0 - (calpha * (mesh_parameters_all_levels[lev])**(-alpha))/tol # use the ComputeTheta function to evaluate theta for level lev
            self.ComputeTheta(lev) # this function bounds theta with a maximum and a minimum value that we set in settings
            theta_i = self.theta_i
            if (theta_i > 0.0) and (theta_i < 1.0):
                coeff2 = np.sum(np.sqrt(np.multiply(model_cost[0:lev+1],BayesianVariance[0:lev+1])))
                coeff2 = coeff2**2.0
                coeff1 = (cphi/(theta_i*tol))**2.0 # formula in case QoI is scalar
            else:
                raise Exception ("The splitting parameter theta_i assumed a value outside the range (0,1) : ",self.theta_i)
            Wtot = coeff1 * coeff2
            # print("print level and correspondent cost",lev,Wtot)
            '''change number of levels if the cost condition is satisfied'''
            if Wtot < Wmin:
                Wmin = Wtot
                current_number_levels = lev
        '''add maximum one level per time'''
        if current_number_levels > Lmin:
            current_number_levels = Lmin + 1
        '''update length of Bayesian variance with respect to the current levels'''
        BayesianVariance = BayesianVariance[0:current_number_levels+1]
        '''update variables'''
        self.BayesianVariance = BayesianVariance
        self.current_number_levels = current_number_levels
        self.previous_number_levels = Lmin
        del(tol,Wmin,current_number_levels,cgamma,gamma,calpha,alpha,cphi,mesh_parameters_all_levels,Lmax,Lmin)

    '''
    function computing the splitting parameter theta \in (0,1)
    '''
    def ComputeTheta(self,level):
        calpha = self.rates_error["calpha"]
        alpha = self.rates_error["alpha"]
        tol = self.tolerance_i
        mesh_param = self.mesh_parameters[level]
        self.theta_i = 1.0 - (calpha * (mesh_param)**(-alpha))/tol
        if (self.theta_i > self.settings["splitting_parameter_max"].GetDouble()):
            self.theta_i = self.settings["splitting_parameter_max"].GetDouble()
        elif (self.theta_i < self.settings["splitting_parameter_min"].GetDouble()):
            self.theta_i = self.settings["splitting_parameter_min"].GetDouble()
        del(calpha,alpha,tol,mesh_param)

    '''
    function computing the updated number of samples for each level for i^th iteration of CMLMC algorithm
    '''
    def ComputeNumberSamples(self):
        current_number_levels = self.current_number_levels
        BayesianVariance = self.BayesianVariance
        min_samples_add = np.multiply(np.ones(current_number_levels+1),self.settings["minimum_add_level"].GetDouble())
        cgamma = self.rates_error["cgamma"]
        gamma  = self.rates_error["gamma"]
        cphi = self.settings["cphi"].GetDouble()
        mesh_parameters_current_levels = self.mesh_parameters[0:current_number_levels+1]
        theta = self.theta_i
        tol = self.tolerance_i
        nsam = self.number_samples
        '''compute optimal number of samples and store previous number of samples'''
        coeff1 = (cphi/(theta*tol))**2.0
        model_cost = np.multiply(cgamma,np.power(mesh_parameters_current_levels,gamma))
        coeff2 = np.sqrt(np.divide(BayesianVariance,model_cost))
        coeff3 = np.sum(np.sqrt(np.multiply(model_cost,BayesianVariance)))
        opt_number_samples = np.multiply(coeff1*coeff3,coeff2)
        # print("optimal number of samples computed = ",opt_number_samples)
        if (len(opt_number_samples) != current_number_levels+1):
            raise Exception ("length optimal number of samples and current optimal number of level not coherent")
        for lev in range (current_number_levels+1):
            opt_number_samples[lev] = np.ceil(opt_number_samples[lev])
            opt_number_samples[lev] = opt_number_samples[lev].astype(int)
        if len(nsam) < len(opt_number_samples):
            for _ in range (len(opt_number_samples)-len(nsam)):
                nsam.append(0)
        '''copy local string and avoid reference since we modify it'''
        previous_number_samples = nsam[:]
        '''evaluate difference between new and previous number of samples'''
        diff_nsam = []
        for lev in range(current_number_levels+1):
            diff_nsam.append(opt_number_samples[lev] - nsam[lev])
            '''set here that if the optimal number of samples is smaller than
            the previous number of samples, keep the previous number of samples,
            i.e. do not decrease the number of samples w.r.t. previous iteration'''
            if diff_nsam[lev] <= 0.:
                diff_nsam[lev] = 0.
                opt_number_samples[lev] = nsam[lev]
            '''set here to have an addition of minimum min_samples_add per level'''
            if (diff_nsam[lev] > 0.) and (diff_nsam[lev] < min_samples_add[lev]):
                diff_nsam[lev] = min_samples_add[lev]
                opt_number_samples[lev] = nsam[lev] + diff_nsam[lev]
            nsam[lev] = opt_number_samples[lev]
        '''convert current number of samples and difference number of samples to integers'''
        for lev in range (current_number_levels+1):
            diff_nsam[lev] = int(diff_nsam[lev])
            nsam[lev] = int(nsam[lev])
        '''update variables and delete local variables'''
        self.number_samples = nsam
        self.difference_number_samples = diff_nsam
        self.previous_number_samples = previous_number_samples
        del(current_number_levels,BayesianVariance,min_samples_add,cgamma,gamma,cphi,mesh_parameters_current_levels,\
        theta,tol,nsam,opt_number_samples,previous_number_samples,diff_nsam)

    '''
    function computing the mlmc estimator for the mean of the Quantity of Interest
    '''
    def ComputeMeanMLMCQoI(self):
        self.mean_mlmc_QoI = np.sum(self.difference_QoI.mean)

    '''
    function computing the total error:
    TErr = bias contrinution + statistical error contribution
    bias contribution B ~= abs(E^MC[QoI_{L}-QoI_{L-1}])
    statistical error contribution SE = \sum_{i=0}^{L}(Var^MC[Y_l]/N_l)
    '''
    def ComputeTotalErrorMLMC(self):
        self.difference_QoI.bias_error = np.abs(self.difference_QoI.mean[self.current_number_levels])
        variance_from_bayesian = np.zeros(np.size(self.number_samples))
        for lev in range(self.current_number_levels+1):
            variance_from_bayesian[lev] = self.BayesianVariance[lev]/self.number_samples[lev]
        self.difference_QoI.statistical_error = self.settings["cphi"].GetDouble() * np.sqrt(np.sum(variance_from_bayesian))
        TErr = self.difference_QoI.bias_error + self.difference_QoI.statistical_error
        self.TErr = TErr


class MultilevelMonteCarloResults(object):
    '''The base class for the MultilevelMonteCarloResults-classes'''
    def __init__(self):
        '''constructor of the MultilevelMonteCarloResults-Object
        Keyword arguments:
        self : an instance of a class
        '''
        '''Quantity of Interest'''
        self.QoI = []
        '''time cost'''
        self.time_ML = []
        '''level of QoI and time_ML'''
        self.finer_level = 0