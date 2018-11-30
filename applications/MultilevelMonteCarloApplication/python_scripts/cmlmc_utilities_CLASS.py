from __future__ import absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import numpy as np
import KratosMultiphysics
import time

# Import exaqute
# from exaqute.ExaquteTaskPyCOMPSs import *   # to exequte with pycompss
# from exaqute.ExaquteTaskHyperLoom import *  # to exequte with the IT4 scheduler
# from exaqute.ExaquteTaskLocal import *      # to execute with python3
# get_value_from_remote is the equivalent of compss_wait_on
# in the future, when everything is integrated with the it4i team, putting exaqute.ExaquteTaskHyperLoom you can launch your code with their scheduler instead of BSC


'''
This utility contains all the functions to perform the Continuation Multilevel Monte Carlo (CMLMC) algorithm described in [PNL17]

References:
[PNL17] M. Pisaroni; F. Nobile; P. Leyland : A Continuation Multi Level Monte Carlo (C-MLMC) method for uncertainty quantification in compressible inviscid aerodynamics; Computer Methods in Applied Mechanics and Engineering, vol 326, pp 20-50, 2017. DOI : 10.1016/j.cma.2017.07.030.
'''

class StatisticalVariable(object):
    '''The base class for the quantity of interest and other statistical variables computed'''
    def __init__(self):
        '''
        Constructor of the class

        Keyword arguments:
        self : an instance of a class
        '''

        '''values of the variable, divided per level'''
        self.values = []
        '''mean of the variable per each level'''
        self.mean = []
        '''sample variance of the variable per each level'''
        self.sample_variance = []
        '''second moment of the variable per each level'''
        self.second_moment = []
        '''bias error of the variable'''
        self.bias_error = None
        '''statistical error of the variable'''
        self.statistical_error = None

    '''
    update mean and second moment values
    M_{2,n} = sum_{i=1}^{n} (x_i - mean(x)_n)^2
    M_{2,n} = M_{2,n-1} + (x_n - mean(x)_{n-1}) * (x_n - mean(x)_{n})
    s_n^2 = M_{2,n} / (n-1)
    '''
    def UpdateOnepassMeanVariance(self,level,i_sample):
        sample = self.values[level][i_sample]
        old_mean = self.mean[level]
        old_M2 = self.second_moment[level]
        nsamples = i_sample + 1
        delta = np.subtract(sample, old_mean)
        if nsamples == 1:
            new_mean = sample
            new_M2 = np.zeros(np.size(sample))
            new_M2 = np.asscalar(new_M2) # do so to have a list of scalars, and not a list of arrays of one element
            new_sample_variance = np.zeros(np.size(sample))
            new_sample_variance = np.asscalar(new_sample_variance) # do so to have a list of scalars, and not a list of arrays of one element
        else:
            new_mean = old_mean + np.divide(delta,nsamples)
            new_M2 = old_M2 + delta*np.subtract(sample,new_mean)
            new_sample_variance = np.divide(new_M2,np.subtract(nsamples,1))
        self.mean[level] = new_mean
        self.second_moment[level] = new_M2
        self.sample_variance[level] = new_sample_variance
        del(new_mean, new_M2, new_sample_variance)


class MultilevelMonteCarlo(object):
    '''The base class for the MultilevelMonteCarlo-classes'''
    def __init__(self,settings):
        '''The constructor of the MultilevelMonteCarlo-Object

        Keyword arguments:
        self     : an instance of a class
        settings : the settings of the Multilevel Monte Carlo simulation
        '''

        '''
        k0   : Certainty Parameter 0 rates
        k1   : Certainty Parameter 1 rates
        r1   : Cost increase first iterations C-MLMC
        r2   : Cost increase final iterations C-MLMC
        tol0 : Tolerance iter 0
        tolF : Tolerance final
        cphi : Confidence on tolerance
        N0   : Number of samples for iter 0
        L0   : Number of levels for iter 0
        Lmax : Maximum number of levels
        mesh_refinement_coefficient : coefficient of mesh refinement
        initial_mesh_size : size of first mesh considered
        minimum_add_level : minimum number of samples to add if at least one is going to be added
        '''
        self.settings = {"k0":settings[0],\
                         "k1":settings[1],\
                         "r1":settings[2],\
                         "r2":settings[3],\
                         "tol0":settings[4],\
                         "tolF":settings[5],\
                         "cphi":settings[6],\
                         "N0":settings[7],\
                         "Lscreening":settings[8],\
                         "Lmax":settings[9],\
                         "mesh_refinement_coefficient":settings[10],\
                         "initial_mesh_size":settings[11]}
        self.settings["minimum_add_level"] = 6.
        '''current_number_levels : number of levels of current iteration'''
        self.current_number_levels = self.settings["Lscreening"]
        '''previous_number_levels : number of levels of previous iteration'''
        self.previous_number_levels = None
        '''number_samples : total number of samples at current iteration'''        
        self.number_samples = [self.settings["N0"] for i in range (self.settings["Lscreening"]+1)]
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
        self.size_mesh = []
        '''BayesianVariance : Bayesian variance'''
        self.BayesianVariance = []
        '''number_iterations : theoretical number of iterations the MLMC algorithm will perform'''
        self.number_iterations_iE = None
        '''convergence : boolean variable defining if MLMC algorithm is convergenced'''
        self.convergence = False
        '''current_iteration : current iteration of MLMC algorithm'''
        self.current_iteration = None
        '''tolerance_i : tolerance of i^th-iteration considered in MLMC algorithm'''
        self.tolerance_i = None
        '''theta_i : splitting parameter \in(0,1) that affects bias and statistical error in the computation of the total error'''
        self.theta_i = None
        '''mean_mlmc_QoI : MLMC estimator for the mean value of the Quantity of Interest'''
        self.mean_mlmc_QoI = None
        '''TErr : total error of MLMC algorithm, the sum of bias and statistical error is an overstmation of the real total error
                  TErr := \abs(E^MLMC[QoI] - E[QoI])'''
        self.TErr = None

        '''difference_QoI : Quantity of Interest of the considered problem organized in consecutive levels
                            difference_QoI.values := Y_l = QoI_M_l - Q_M_l-1'''
        self.difference_QoI = StatisticalVariable()
        self.difference_QoI.values = [[] for i in range (self.settings["Lscreening"]+1)] # list containing Y_{l}^{i} = Q_{m_l} - Q_{m_{l-1}}
        '''time_ML : time to perform a single MLMC simulation (i.e. one value of difference_QoI.values) organized in consecutive levels'''
        self.time_ML = StatisticalVariable()
        self.time_ML.values = [[] for i in range (self.settings["Lscreening"]+1)] # list containing the time to compute the level=l simulations

        '''########################################################################
        # observation: levels start from level 0                                  #
        #              length arrays and lists starts from 1                      #
        # then we have a difference of 1 between length lists and levels          #
        # e.g. self.current_level = len(self.number_samples) - 1                  #
        #      or                                                                 #
        #      self.current_level = len(self.difference_QoI.difference_value) - 1 #
        ########################################################################'''


    '''
    function finalizing the screening phase of the MLMC algorithm
    Usage: It is designed to be called ONCE, AFTER the screening phase
    '''
    def FinalizeScreeningPhase(self):
        '''prepare lists'''
        self.difference_QoI.mean = [[] for i in range (self.settings["Lscreening"]+1)]
        self.difference_QoI.sample_variance = [[] for i in range (self.settings["Lscreening"]+1)]
        self.difference_QoI.second_moment = [[] for i in range (self.settings["Lscreening"]+1)]
        self.time_ML.mean = [[] for i in range (self.settings["Lscreening"]+1)]
        self.time_ML.sample_variance = [[] for i in range (self.settings["Lscreening"]+1)]
        self.time_ML.second_moment = [[] for i in range (self.settings["Lscreening"]+1)]
        '''compute mean, sample variance and second moment for difference QoI and time ML'''
        for level in range (self.current_number_levels+1):
            for i_sample in range(self.number_samples[level]):
                self.difference_QoI.UpdateOnepassMeanVariance(level,i_sample)
                self.time_ML.UpdateOnepassMeanVariance(level,i_sample)
        '''compute mesh parameter for each mesh'''
        self.ComputeMeshParameters()
        '''compute parameters by least square fit to estimate Bayesian VAR'''
        self.ComputeRatesLS()
        '''compute Bayesian VAR V^c[Y_l]'''
        self.EstimateBayesianVariance(self.current_number_levels)
        '''compute i_E, number of iterations'''
        self.ComputeNumberIterationsMLMC()
        '''start first iteration, we are entering in the MLMC algorithm'''
        self.current_iteration = 1
    
    '''
    function performing all the required operations that should be executed
    (for each step) BEFORE the MLMC solution step
    '''
    def InitializeMLMCPhase(self):
        '''compute tolerance for the i^th iteration'''
        self.ComputeTolerancei()
        '''Compute Optimal Number of Levels for iteration i L_i'''
        self.ComputeLevels()
        '''compute theta splitting parameter according to the current_number_levels and tolerance_i'''
        self.ComputeTheta()
        '''compute number of samples according to BayesianVariance and theta_i parameters'''
        self.ComputeNumberSamples()
        '''prepare lists'''
        for i in range (self.current_number_levels - self.previous_number_levels): # append a list for the new level
            self.difference_QoI.values.append([])
            self.time_ML.values.append([])


    '''function performing all the required operations that should be executed
    (for each step) AFTER the MLMC solution step'''
    def FinalizeMLMCPhase(self):
        '''prepare lists'''
        for i in range (self.current_number_levels - self.previous_number_levels): # append a list for the new level
            self.difference_QoI.mean.append([])
            self.difference_QoI.sample_variance.append([])
            self.difference_QoI.second_moment.append([])
            self.time_ML.mean.append([])
            self.time_ML.sample_variance.append([])
            self.time_ML.second_moment.append([])
        '''compute mean, second moment and sample variance'''
        for level in range (self.current_number_levels+1):
            # for i_sample in range(self.difference_number_samples[level]):
            for i_sample in range(self.previous_number_samples[level],self.number_samples[level]):
                self.difference_QoI.UpdateOnepassMeanVariance(level,i_sample)
                self.time_ML.UpdateOnepassMeanVariance(level,i_sample)
        '''compute estimatior MLMC mean QoI'''
        self.compute_mean_mlmc_QoI()
        '''compute parameters by least square fit'''
        self.ComputeRatesLS()
        '''compute Bayesian variance'''
        self.EstimateBayesianVariance(self.current_number_levels)
        '''compute total error of the MLMC simulation'''
        self.ComputeTotalErrorMLMC()
        '''update number of levels'''
        self.previous_number_levels = self.current_number_levels
        '''convergence reached if: i) current_iteration >= number_iterations_iE
                                  ii) TErr < tolerance_i
           if not update current_iteration'''
        if (self.current_iteration >= self.number_iterations_iE) and (self.TErr < self.tolerance_i):
            self.convergence = True
        else:
            self.current_iteration = self.current_iteration + 1


    '''
    function giving as output the mesh discretization parameter
    the mesh parameter is the reciprocal of the minimum mesh size of the grid
    h_lev=h_0*M^(-lev)
    '''
    def ComputeMeshParameters(self):
        h0 = self.settings["initial_mesh_size"]
        M  = self.settings["mesh_refinement_coefficient"]
        for level in range(self.settings["Lmax"]+1):
            h_current_level = h0 * M**(-level)
            mesh_parameter_current_level = h_current_level**(-1)
            self.size_mesh.append(h_current_level)
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
        '''update the dictionary'''
        self.rates_error["calpha"] = C1
        self.rates_error["alpha"] = alpha
        self.rates_error["cbeta"] = C2
        self.rates_error["beta"] = beta
        self.rates_error["cgamma"] = C3
        self.rates_error["gamma"] = gamma
        del(bias_ratesLS,variance_ratesLS,cost_ML_ratesLS,mesh_param_ratesLS,C1,C2,C3,alpha,beta,gamma)


    '''
    function performing the Bayesian update of the variance
    using samples generated on all levels in order to locally improve the estimation of Var[difference_QoI]
    '''
    def EstimateBayesianVariance(self,levels): # need to keep levels because in ComputeLevels I use the maximum number of levels
        '''use local variables'''
        k0 = self.settings["k0"]
        k1 = self.settings["k1"]
        calfa = self.rates_error["calpha"]
        alfa  = self.rates_error["alpha"]
        cbeta = self.rates_error["cbeta"]
        beta  = self.rates_error["beta"]
        mesh_param = self.mesh_parameters
        '''use local variables, in order to not modify the global variables'''
        mean_local = self.difference_QoI.mean[:]
        variance_local = self.difference_QoI.sample_variance[:]
        nsam_local = self.number_samples[:]

        if len(mean_local) < (levels+1):
            for i in range (0,(levels+1)-len(mean_local)):
                mean_local.append(0.0)
        if len(variance_local) < (levels+1):
            for i in range (0,(levels+1)-len(variance_local)):
                variance_local.append(0.0)
        if len(nsam_local) < (levels+1):
            for i in range ((levels+1)-len(nsam_local)):
                nsam_local.append(0)

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
    function computing the iteration parameter i_E
    the first i_E iterations are needed to obtain increasingly accurate estimates of the problem dependent parameters P=[calpha,alpha,cbeta,beta,cgamma,gamma]
    the iterations i>i_E prevent redundant computations due to fluctuations in the estimate of P=[calpha,alpha,cbeta,beta,cgamma,gamma] by solving the problem for a slightly smaller tolerance than the desired one
    the function "compute_tolerance_i" computes the tolerance for itaeration "i"
    '''
    def ComputeNumberIterationsMLMC(self):
        tolF = self.settings["tolF"]
        tol0 = self.settings["tol0"]
        r2 = self.settings["r2"]
        r1 = self.settings["r1"]
        self.number_iterations_iE = np.floor((-np.log(tolF)+np.log(r2)+np.log(tol0))/(np.log(r1)))

    def ComputeTolerancei(self):
        tolF = self.settings["tolF"]
        r2 = self.settings["r2"]
        r1 = self.settings["r1"]
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
    function computing the number of levels for iteration "i" of the cmlmc algorithm
    '''
    def ComputeLevels(self):
        '''observe I have already computed mesh parameters for all the possible levels, i.e. up to Lmax'''
        tol = self.tolerance_i
        Wmin   = 1e10
        Lopt_local = self.current_number_levels
        nsam_local = self.number_samples
        Cgamma = self.rates_error["cgamma"]
        gamma  = self.rates_error["gamma"]
        Calpha = self.rates_error["calpha"]
        alpha = self.rates_error["alpha"]
        Cphi = self.settings["cphi"]
        mesh_param_all = self.mesh_parameters
        mean = self.difference_QoI.mean
        variance = self.difference_QoI.sample_variance
        Lmax = self.settings["Lmax"]
        Lmin = self.current_number_levels

        if len(self.BayesianVariance) < (Lmax+1):
            self.EstimateBayesianVariance(Lmax)
        BayesianVariance = self.BayesianVariance
        '''now both mesh_param_all and BayesianVariance have length = Lmax + 1'''
        model_cost = np.multiply(Cgamma,np.power(mesh_param_all,gamma))
        '''also model_cost has length = Lmax + 1'''
                                                                
        
        for lev in range(Lmin, Lmax+1):
            '''it is not mandatory to increase the number of levels and we may continue using the number of levels of the previous iteration, i.e. Lmin
            consider theta_i = 1.0 - (calpha*(M_L**(-alpha)))/tol_i'''
            theta_i = 1.0 - (Calpha * (mesh_param_all[lev])**(-alpha))/tol # I do not call the def "theta_model",
                                                                    # because I cannot use a task inside a task
            if (theta_i > 0.0) and (theta_i < 1.0):
                '''update the cost in future: for the levels we know use cost_ML
                                            for the levels we do not know use "cgamma*Ml**gamma", where Ml is the mesh parameter'''
                coeff2 = np.sum(np.sqrt(np.multiply(model_cost[0:lev+1],BayesianVariance[0:lev+1])))
                coeff2 = coeff2**2.0
                coeff1 = (Cphi/(theta_i*tol))**2.0 # formula in case QoI is scalar, if QoI use the formula described in [PNL16]
                
            else:
                raise Exception ("The splitting parameter theta_i assumed a value outside the range (0,1)")

            Wtot = coeff1 * coeff2
            print("print level and correspondent cost",lev,Wtot)
            if Wtot < Wmin:
                Wmin = Wtot
                Lopt_local = lev
        
        if Lopt_local > Lmin:
            Lopt_local = Lmin + 1
            '''i.e. add one level per time
            note that the number of levels starts from 0, not from 1,
            so we have a difference of one between the number of levels
            and the length of the arrays
            (e.g. difference_value, mesh_param_all or number_sample)'''
        BayesianVariance = BayesianVariance[0:Lopt_local+1]
        '''need to leave Lopt, and so the new BayesianVariance value,
        because I need this value to compute number of samples'''

        self.BayesianVariance = BayesianVariance
        self.current_number_levels = Lopt_local
        self.previous_number_levels = Lmin


    '''
    function computing the splitting parameter theta \in (0,1)
    '''
    def ComputeTheta(self):
        Calpha = self.rates_error["calpha"]
        alpha = self.rates_error["alpha"]
        tol = self.tolerance_i
        mesh_param = self.mesh_parameters[self.current_number_levels]
        self.theta_i = 1.0 - (Calpha * (mesh_param)**(-alpha))/tol
        if (self.theta_i < 0.0) or (self.theta_i > 1.0):
            raise Exception ("The splitting parameter theta_i assumed a value outside the range (0,1)")

    '''
    function computing the new number of samples for each level for the iteration "i" of the cmlmc algorithm
    '''
    def ComputeNumberSamples(self):
        L_opt = self.current_number_levels
        BayesianVariance = self.BayesianVariance
        minNadd = np.multiply(np.ones(L_opt+1),self.settings["minimum_add_level"])
        Cgamma = self.rates_error["cgamma"]
        gamma  = self.rates_error["gamma"]
        Cphi = self.settings["cphi"]
        mesh_param_local = self.mesh_parameters[0:L_opt+1]
        theta = self.theta_i
        tol = self.tolerance_i
        nsam = self.number_samples

        coeff1 = (Cphi/(theta*tol))**2.0
        model_cost = np.multiply(Cgamma,np.power(mesh_param_local,gamma))

        coeff2 = np.sqrt(np.divide(BayesianVariance,model_cost))
        coeff3 = np.sum(np.sqrt(np.multiply(model_cost,BayesianVariance)))
        
        opt_number_samples = np.multiply(coeff1*coeff3,coeff2)
        print("optimal number of samples computed = ",opt_number_samples)
        
        for i in range (0,len(opt_number_samples)):
            opt_number_samples[i] = np.ceil(opt_number_samples[i])
            opt_number_samples[i] = opt_number_samples[i].astype(int)

        if len(nsam) < len(opt_number_samples):
            for i in range (0,len(opt_number_samples)-len(nsam)):
                nsam.append(0)
        previous_number_samples = nsam[:]

        dNsam = []
        for l in range(0,L_opt+1):
            dNsam.append(opt_number_samples[l] - nsam[l])
            if dNsam[l] <= 0.:
                dNsam[l] = 0.
                opt_number_samples[l] = nsam[l]
                '''i.e. here I set that if NlOPT[l] is smaller than the previous
                number of samples, I keep the previous number of samples'''
    
            if (dNsam[l] > 0.) and (dNsam[l] < minNadd[l]):
                dNsam[l] = minNadd[l]
                opt_number_samples[l] = nsam[l] + dNsam[l]
                '''i.e. the minimum addition of samples is given by the array minNadd
                so if the new number of samples would be smaller than the old
                number of samples, I set to have an addition of minNadd[l] samples'''
            
            nsam[l] = opt_number_samples[l]
        for i in range (0,len(dNsam)):
            dNsam[i] = int(dNsam[i])
            # dNsam[i] = dNsam[i].astype(int)
            nsam[i] = int(nsam[i])
            # nsam[i] = nsam[i].astype(int)

        '''note that I do not decrease the number of samples wrt previous MLMC iterations!!'''
        self.number_samples = nsam
        self.difference_number_samples = dNsam
        self.previous_number_samples = previous_number_samples
        return nsam,dNsam,previous_number_samples

    '''
    function computing the mlmc estimator for the mean of the Quantity of Interest
    '''
    def compute_mean_mlmc_QoI(self):
        self.mean_mlmc_QoI = np.sum(self.difference_QoI.mean)


    '''
    function computing the total error:
    TErr = bias contrinution + statistical error contribution
    bias contribution B ~= abs(E^MC[Q_{L}-Q_{L-1}])
    statistical error contribution SE = \sum_{i=0}^{L}(Var^MC[Y_l]/N_l)
    '''
    def ComputeTotalErrorMLMC(self):
        self.difference_QoI.bias_error = np.abs(self.difference_QoI.mean[self.current_number_levels])
        var_bayes = np.zeros(np.size(self.number_samples))
        for i in range(self.current_number_levels+1):
            var_bayes[i] = self.BayesianVariance[i]/self.number_samples[i]
        self.difference_QoI.statistical_error = self.settings["cphi"] * np.sqrt(np.sum(var_bayes))
        TErr = self.difference_QoI.bias_error + self.difference_QoI.statistical_error
        self.TErr = TErr
