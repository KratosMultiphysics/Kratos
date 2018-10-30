from __future__ import absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import numpy as np
# import time

'''
This file contains all the functions to perform the Continuation Multilevel Monte Carlo (CMLMC) algorithm described in [PNL17]

References:
[PNL17] M. Pisaroni; F. Nobile; P. Leyland : A Continuation Multi Level Monte Carlo (C-MLMC) method for uncertainty quantification in compressible inviscid aerodynamics; Computer Methods in Applied Mechanics and Engineering, vol 326, pp 20-50, 2017. DOI : 10.1016/j.cma.2017.07.030.
'''


'''
function that gives as output a list containing the mesh discretization parameter
we consider the number of elements of a uniform grid on the same domain to be this parameter
refinement strategy:
uniform mesh on level "lev" with h_lev=(1/N0)*2^(-lev)
i.e. level = 0, h_{lev=0} = 0.25
        level = 1, h_{lev=1} = 0.2/2 = 0.125
        level = 2, h_{lev=2} = 0.2/4 = 0.0625
        ...
'''
def Nf_law(lev):
    # I wrote the following:
    N0 = 4.
    M  = 2.
    NFF = (N0*np.power(M,lev))
    Nf2 = 2*NFF**2
    # Nobile wrote the following:
    # N0 = 5.
    # M  = 2.
    # NFF = (N0*np.power(M,lev))
    # Nf2 = NFF**2
    '''
    NFF is the number of elements on a boundary line
    Nf2 is the number of triangular elements in the square domain (approximately, if mesh non uniform)
    '''
    return Nf2


'''
update mean and second moment values
M_{2,n} = sum_{i=1}^{n} (x_i - mean(x)_n)^2
M_{2,n} = M_{2,n-1} + (x_n - mean(x)_{n-1}) * (x_n - mean(x)_{n})
s_n^2 = M_{2,n} / (n-1)
'''
def update_onepass_M(sample, old_mean, old_M2, nsam):
    delta = np.subtract(sample, old_mean)
    if nsam == 1:
        new_mean = sample
        new_M2 = np.zeros(np.size(sample))
        new_M2 = np.asscalar(new_M2)
        new_sample_variance = np.zeros(np.size(sample))
        new_sample_variance = np.asscalar(new_sample_variance)
        '''do so to have a list of scalars, and not a list of arrays of one element'''
    else:
        new_mean = old_mean + np.divide(delta,nsam)
        new_M2 = old_M2 + delta*np.subtract(sample,new_mean)
        new_sample_variance = compute_sample_variance_from_M2(new_M2,nsam)
    return new_mean, new_M2, new_sample_variance

def compute_sample_variance_from_M2(M2,nsam):
    sample_variance = np.divide(M2,np.subtract(nsam,1))
    return sample_variance


'''
function performing the Bayesian update of the variance using samples generated on all levels in order to locally improve the estimation of Var[Y_l]
see [PNL17] pp. 8-9
'''
def EstimateBayesianVariance(mean,variance,settings_ML,ratesLS,nDoF,nsam,level_local):
    k0 = settings_ML[0]
    k1 = settings_ML[1]
    Calfa = ratesLS[0]
    alfa  = ratesLS[1]
    Cbeta = ratesLS[2]
    beta  = ratesLS[3]

    '''use local variables, in order to not modify the global variables'''
    mean_local = mean[:]
    variance_local = variance[:]
    nsam_local = nsam[:]
    if len(mean_local) < (level_local+1):
        for i in range (0,(level_local+1)-len(mean_local)):
            mean_local.append(0.0)
    if len(variance_local) < (level_local+1):
        for i in range (0,(level_local+1)-len(variance_local)):
            variance_local.append(0.0)
    if len(nsam_local) < (level_local+1):
        for i in range (0,(level_local+1)-len(nsam_local)):
            nsam_local.append(0)

    BayesianVariance = []
    for level in range (0, (level_local+1)):
        mu = Calfa*nDoF[level]**(-alfa)
        lam = (1/Cbeta)*nDoF[level]**(beta)
        G1_l = 0.5 + np.multiply(k1,lam) + np.divide(nsam_local[level],2.0)
        G2_l = k1 + (nsam_local[level]-1)*0.5*variance_local[level] + k0*nsam_local[level]*((mean_local[level]-mu)**2)/(2.0*(k0+nsam_local[level]))
        BayesianVariance.append(np.divide(G2_l,G1_l-0.5))
    return BayesianVariance


'''
function computing the iteration parameter i_E
the first i_E iterations are needed to obtain increasingly accurate estimates of the problem dependent parameters P=[calpha,alpha,cbeta,beta,cgamma,gamma]
the iterations i>i_E prevent redundant computations due to fluctuations in the estimate of P=[calpha,alpha,cbeta,beta,cgamma,gamma] by solving the problem for a slightly smaller tolerance than the desired one
the function "compute_tolerance_i" computes the tolerance for itaeration "i"
see [PNL16] pp. 7-8
'''
def compute_iE_cmlmc(settings_ML):
    tolF = settings_ML[5]
    tol0 = settings_ML[4]
    r2 = settings_ML[3]
    r1 = settings_ML[2]
    iE_cmlmc = np.floor((-np.log(tolF)+np.log(r2)+np.log(tol0))/(np.log(r1)))
    return iE_cmlmc

def compute_tolerance_i(settings_ML,iE,iter_def):
    r1 = settings_ML[2]
    r2 = settings_ML[3]
    tolF = settings_ML[5]
    if iter_def <= iE:
        tol = (r1**(iE-iter_def) * r2**(-1))*tolF
    else:
        tol = (r2**(iE-iter_def) * r2**(-1))*tolF
    return tol


'''
function computing the problem dependent parameters P=[calpha,alpha,cbeta,beta,cgamma,gamma] using least squares fit
see [PNL17] pp.7-8
'''
def compute_ratesLS(bias_ratesLS,variance_ratesLS,cost_ML_ratesLS,ndof_ratesLS):
    bias_ratesLS = np.abs(bias_ratesLS)

    '''##################### MEAN - alpha ########################################
    NUMPY: linear fit
    why not considered also M_{L=0}?'''
    pa = np.polyfit(np.log2(ndof_ratesLS[1::]),np.log2(bias_ratesLS[1::]),1) # in [PNL17] not used level = 0
    alpha   = -pa[0]
    C1      = 2**pa[1]

    '''##################### VAR - beta ##########################################
    NUMPY: linear fit
    why not considered also M_{L=0}?'''
    pb          = np.polyfit(np.log2(ndof_ratesLS[1::]),np.log2(variance_ratesLS[1::]),1) # in [PNL17] not used level = 0
    beta        = -pb[0]
    C2          = 2**pb[1]

    '''################# COST - gamma ############################################ 
    NUMPY: linear fit'''
    pg          = np.polyfit(np.log2(ndof_ratesLS),np.log2(cost_ML_ratesLS),1)
    gamma       = pg[0]
    C3          = 2**pg[1]

    paramLS=[C1,alpha,C2,beta,C3,gamma]
    return paramLS


'''
function computing the splitting parameter theta \in (0,1)
see [PNL17] pp. 5-6,8-11
'''
def theta_model(ratesLS,toll,nDoF):
    Calpha = ratesLS[0]
    alpha = ratesLS[1]
    theta_i = 1.0 - (Calpha * (nDoF)**(-alpha))/toll
    return theta_i


'''
function computing the number of levels for iteration "i" of the cmlmc algorithm
see [PNL17] pp. 8-11
'''
def compute_levels(tol,nsam,ratesLS,ndof_all,BayesianVariance,mean,variance,settings_ML,Lmax,Lmin):
    '''observe I have already computed ndof_all for all the possible levels, i.e. up to Lmax'''
    Wmin   = 1e10
    Lopt_local = Lmin
    nsam_local = nsam
    Cgamma = ratesLS[4]
    gamma  = ratesLS[5]
    Calpha = ratesLS[0]
    alpha = ratesLS[1]
    Cphi = settings_ML[6]

    if len(BayesianVariance) < (Lmax+1):
        BayesianVariance = EstimateBayesianVariance(mean,variance,settings_ML,ratesLS,ndof_all,nsam_local,Lmax)
    '''now both ndof_all and BayesianVariance have length = Lmax + 1'''
    model_cost = np.multiply(Cgamma,np.power(ndof_all,gamma))
    '''also model_cost has length = Lmax + 1'''
                                                              
    
    for lev in range(Lmin, Lmax+1):
        '''it is not mandatory to increase the number of levels and we may continue using the number of levels of the previous iteration, i.e. Lmin
        consider theta_i = 1.0 - (calpha*(M_L**(-alpha)))/tol_i'''
        theta_i = 1.0 - (Calpha * (ndof_all[lev])**(-alpha))/tol # I do not call the def "theta_model",
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
        # print("print level and correspondent cost",lev,Wtot)
        if Wtot < Wmin:
            Wmin = Wtot
            Lopt_local = lev
    
    if Lopt_local > Lmin:
        Lopt_local = Lmin + 1
        '''i.e. add one level per time
        note that the number of levels starts from 0, not from 1,
        so we have a difference of one between the number of levels
        and the length of the arrays
        (e.g. difference_value, ndof_all or number_sample)'''
    BayesianVariance = BayesianVariance[0:Lopt_local+1]
    '''need to leave Lopt, and so the new BayesianVariance value,
    because I need this value to compute number of samples'''
    
    return Lopt_local, BayesianVariance, Lmin


'''
function computing the new number of samples for each level for the iteration "i" of the cmlmc algorithm
see [PNL17] pp. 6,10
'''
def compute_number_samples(L_opt,BayesianVariance,ratesLS,theta,tol,nDoF,nsam,settings_ML):
    minNadd = np.multiply(np.ones(L_opt+1),6.)
    Cgamma = ratesLS[4]
    gamma  = ratesLS[5]
    Cphi = settings_ML[6]
    ndof_local = nDoF[0:L_opt+1]

    coeff1 = (Cphi/(theta*tol))**2.0
    model_cost = np.multiply(Cgamma,np.power(ndof_local,gamma))

    coeff2 = np.sqrt(np.divide(BayesianVariance,model_cost))
    coeff3 = np.sum(np.sqrt(np.multiply(model_cost,BayesianVariance)))
       
    opt_number_samples = np.multiply(coeff1*coeff3,coeff2)
    
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

    print("new number of samples = ",nsam,"difference with previous iteration = ",dNsam,"old number samples = ",previous_number_samples)
    '''note that I do not decrease the number of samples wrt previous MLMC iterations!!'''
    return nsam,dNsam,previous_number_samples


'''
function computing the mlmc estimator for the mean of the Quantity of Interest
see [PNL17] p. 4
'''
def compute_mean_mlmc_QoI(mean_array):
    mean_mlmc = np.sum(mean_array)
    return mean_mlmc


'''
function computing the total error:
TErr = bias contrinution + statistical error contribution
bias contribution B ~= abs(E^MC[Q_{L}-Q_{L-1}])
statistical error contribution SE = \sum_{i=0}^{L}(Var^MC[Y_l]/N_l)
see [PNL17] pp. 3-7
'''
def compute_total_error_MLMC(mean_difference_QoI,number_samples,L_opt,BayesianVariance,settings_ML):
    bias_error = np.abs(mean_difference_QoI[L_opt])
    var_bayes = np.zeros(np.size(number_samples))
    for i in range(0,L_opt+1):
        var_bayes[i] = BayesianVariance[i]/number_samples[i]
    TErr = bias_error + settings_ML[6]*np.sqrt(np.sum(var_bayes))
    return TErr