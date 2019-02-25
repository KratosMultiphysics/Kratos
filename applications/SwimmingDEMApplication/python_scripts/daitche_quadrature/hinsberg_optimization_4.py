import math 
import cmath
import mpmath
import matplotlib.pyplot as plt
from bigfloat import *
import numpy as np
import scipy
import scipy.optimize as opt
import scipy.integrate as integ
import cmath

def DawsonIntegral(x):
    return 0.5 * sqrt(math.pi) * exp(- x ** 2) * scipy.special.erfi(x)

def FindZero(f, x0):
    return opt.newton(f, 0.5)

class K_B:
    def __init__(self):
        pass
    def f(self, t):
        return t ** (- 0.5)
    
    def df(self, t):
        return - 0.5 * t ** (- 1.5)

class K_component:
    def __init__(self, taus, ai, i):
        self.taus = taus
        self.position = i
        self.t = sum([math.exp(tau) for tau in taus[:i + 1]])
        self.sqrt_e = math.sqrt(math.exp(1))
        self.alpha = self.sqrt_e * self.t ** (- 0.5)
        self.beta = - 0.5 / self.t
        self.a = ai
        
    def alpha_prime(self, i):
        k = self.position
        
        if i > k:
            return 0.
        else:
            return - 0.5 * self.sqrt_e * self.t ** (- 1.5) * math.exp(self.taus[i])
        
    def alpha_prime_2(self, i, j):
        k = self.position
        if i > k or j > k:
            return 0.
        else:
            d2dt = 3. / 4 * self.sqrt_e * self.t ** (- 2.5) * math.exp(self.taus[i] + self.taus[j])
            if i == j:
                d2dt += - 0.5 * self.sqrt_e * self.t ** (- 1.5) * math.exp(self.taus[j])
            return d2dt
        
    def beta_prime(self, i):
        k = self.position
        if i > k:
            return 0.
        else:
            return 0.5 * self.t ** (- 2.) * math.exp(self.taus[i])
      
    def beta_prime_2(self, i, j):
        k = self.position
        if i > k or j > k:
            return 0.
        else:
            d2dt = - self.t ** (- 3.) * math.exp(self.taus[i] + self.taus[j])
            if i == j:
                d2dt += 0.5 * self.t ** (- 2.) * math.exp(self.taus[j])
            return d2dt
      
    def f(self, t):
        a = self.alpha
        b = self.beta
        return a * exp(b * t)
    
    def df(self, t):
        a = self.alpha
        b = self.beta
        return b * a * exp(b * t)

class K_sum:
    def __init__(self, ais, tis):
        self.tis = tis
        self.Ks = [K_component(tis, ais[i], i) for i in range(len(tis))]
        self.m = len(tis)
        
    def f(self, t):    
        Ks = self.Ks
        tis = self.tis
        m = self.m 
        return sum([Ks[i].a * Ks[i].f(t) for i in range(m)])
    
    def df(self, t):
        Ks = self.Ks
        tis = self.tis
        m = self.m 
        return sum([Ks[i].a * Ks[i].df(t) for i in range(m)])

def ApproximateQuadrature(times, f):
    values = [0.0 for t in times]
    acc_sum = 2 * math.sqrt(times[-1] - times[-2]) * f(times[-1])

    for i in range(len(values) - 1):
        if i == 0:
            delta_t = times[1] - times[0]
        else:
            delta_t = times[i] - times[i - 1]
                    
        acc_sum += 0.5 * delta_t * (f(times[i]) + f(times[i - 1]))
        
    return acc_sum

def SubstituteRichardsons(approx_successive_values, k, order, level = - 1):
    with precision(200):
        one = BigFloat(1)
        n = len(approx_successive_values)
        if level > n or level < 0:
            max_n = n
        else:
            max_n = level
        
        richardsons = [value for value in approx_successive_values]
        
        while max_n:
            max_n -= 1
            n -= 1
            new_richardsons = []
            
            for i in range(n):
                new_richardsons.append((k ** order * richardsons[i + 1] - richardsons[i]) / (one * (k ** order - 1)))
                
            richardsons = [value for value in new_richardsons]        
            
            for i in range(n):
                approx_successive_values[- i - 1] = richardsons[- i - 1]
            order += 1

class Functional:
    def __init__(self):
        self.sqrt_pi = sqrt(math.pi)
    
    def Define(self, ais, tis):
        self.K = K_sum(ais, tis)
        self.K_1 = self.K.f(1)
        big_t = 1e3
        self.K_big_t = self.K.f(big_t)
        self.K_B = K_B()
    
    def Fmod(self):
        K = self.K
        K_1 = self.K_1
        K_B = self.K_B
        
        def f(t):            
            return float(t * (- 0.5 / abs(t) ** 1.5 - K.df(t)) ** 2)
        
        fixed_part = (1. - K_1) ** 2
        big_t = 500.
        def Integrate(n):
            times = [1.0 + ((big_t - 1.0) * i) / n for i in range(n + 1)]      
            #print(times)
            
            delta_t = times[1] - times[0]
            acc_sum = 0.0
            for i in range(len(times) - 1):                   
                acc_sum += 0.5 * delta_t * (f(times[i]) + f(times[i - 1]))
            return acc_sum
        early_part_1 = Integrate(16)
        #print("early1", early_part_1)        
        early_part_2 = Integrate(32)
        #print("early2", early_part_2)
        early_part_3 = Integrate(64)
        #print("early3", early_part_3)
        early_part_4 = Integrate(128)
        #print("early4", early_part_4)
        early_part_5 = Integrate(256)
        #print("early5", early_part_5)
        early_part_6 = Integrate(512)   
        #print("early6", early_part_6)
        values = [early_part_1, early_part_2, early_part_3, early_part_4, early_part_5, early_part_6]
        SubstituteRichardsons(values, 2, 2)
        early_part = values[-1]
        return fixed_part + early_part
    
    def F(self):
        K = self.K
        K_1 = self.K_1
        K_B = self.K_B
        
        def f(t):            
            return float(self.K_B.df(t) - K.df(t))
        
        def f_abs(t):
            return abs(f(t))
        
        try:            
            t_cross = FindZero(f, 1.0)
            
            if t_cross <= 1.0:
                print("NO ROOTS BETWEEN 1 AND INFINITY!")
                return  2 * abs(1. - K_1)
            else:
                print("ROOT", t_cross)
        except:
            print("NO ROOTS AT ALL!")
            return 2 * abs(1. - K_1)

        fixed_part = abs(1. - K_1)
        
        K_t_cross = K.f(t_cross)
        early_part = abs(1. / sqrt(t_cross) - 1 + K_1 - K_t_cross)
        late_part = abs(1. / sqrt(t_cross) - K_t_cross)
        residual = fixed_part + early_part + late_part
       
        if residual < 0:
            print("\n NEGATIVE_RES", residual)  
            print("ROOT", t_cross)
            print("smaller t", f(0.9999*t_cross))
            print("bigger t", f(1.00001*t_cross))
            print("standard part", abs(1. - K_1))
            print("bad_pàrt", integral_part)
        return residual
    
        #big_t = 1e1
        #def Integrate(n):
            #times = [1.0 + ((big_t - 1.0) * i) / n for i in range(n + 1)]      
            ##print(times)
            
            #delta_t = times[1] - times[0]
            #acc_sum = 0.0
            #for i in range(len(times) - 1):                   
                #acc_sum += 0.5 * delta_t * (f_abs(times[i]) + f_abs(times[i - 1]))
            #return acc_sum
        #early_part_1 = Integrate(64)
        ##print("early1", early_part_1)
        #early_part_2 = Integrate(128)
        ##print("early2", early_part_2)
        #early_part_3 = Integrate(256)
        ##print("early3", early_part_3)
        #early_part_4 = Integrate(512)   
        ##print("early4", early_part_4)
        #values = [early_part_1, early_part_2, early_part_3, early_part_4]
        #SubstituteRichardsons(values, 2, 2)
        #early_part = values[-1]
        ##print("early", early_part)
        #late_part = 1 / sqrt(big_t) - K.f(big_t)                                

        ##print("fixed", fixed_part)
        ##print("late", late_part)
        #return fixed_part + early_part + late_part
        #para
    


    def dFda(self, i):
        K = self.K
        K_1 = self.K_1
        sqrt_pi = self.sqrt_pi                  
        Ki = K.Ks[i]
        Ki_1 = Ki.f(1)
        alpha_i = Ki.alpha        
        beta_i = Ki.beta        
        sqrt_minus_betai = sqrt(- beta_i)
        sqrt_minus_pi_beta_i = - sqrt_pi * sqrt_minus_betai
        
        non_integral_contribution = - 2 * (1 - K_1) * Ki_1
        Ki_prime_over_sqrt_t = alpha_i * sqrt_minus_pi_beta_i * erfc(sqrt_minus_betai)
        t_Ki_prime_K_prime = 2 * beta_i * Ki_1 * sum([Kk.a * Kk.f(1.) * Kk.beta * (1. - Kk.beta - beta_i) / (Kk.beta + beta_i) ** 2 for Kk in K.Ks])
        
        return  non_integral_contribution + Ki_prime_over_sqrt_t + t_Ki_prime_K_prime    
    
    def dFdt(self, i):        
        K = self.K
        dfdt = 0.
        
        for k in range(i + 1):
            Kk = K.Ks[k]
            alpha_k = Kk.alpha   
            a_k = Kk.a
            dfdt += a_k * (Kk.alpha_prime(i) / alpha_k + Kk.beta_prime(i)) * self.dFda(k)
        return dfdt 
    
    def d2Fda2(self, i, j):
        K = self.K
        Ki = K.Ks[i]
        Kj = K.Ks[j]
        Ki_1 = Ki.f(1)
        Kj_1 = Kj.f(1)
        beta_i = Ki.beta        
        beta_j = Kj.beta        
        gammaij = beta_i * beta_j * (1. - beta_j - beta_i) / (beta_j + beta_i) ** 2 * Ki_1 * Kj_1
        
        return 2 * (Kj_1 * Ki_1 + gammaij)
    
    def d2Fdtda(self, i, j):
        K = self.K
        dtda1 = 0.
        dtda2 = 0.       
        
        for k in range(i + 1):
            Kk = K.Ks[k]
            alpha_k = Kk.alpha   
            a_k = Kk.a
            dtda1 += a_k * (Kk.alpha_prime(i) / alpha_k + Kk.beta_prime(i)) * self.d2Fda2(k, j)
        
        if i >= j:
            Kj = K.Ks[j]
            alpha_j = Kj.alpha   
            dtda2 += (Kj.alpha_prime(i) / alpha_j + Kj.beta_prime(i)) * self.dFda(j)
        
        return dtda1 + dtda2

    def d2Fdt2(self, i, j):
        K = self.K
        dtdt1 = 0.
        dtdt2 = 0.    
        
        for k in range(i + 1):
            Kk = K.Ks[k]
            alpha_k = Kk.alpha   
            a_k = Kk.a
            dtdt1 += a_k * (Kk.alpha_prime(i) / alpha_k + Kk.beta_prime(i)) * self.d2Fdtda(j, k)     

        if i >= j:
            for k in range(i + 1):
                Kk = K.Ks[k]
                alpha_k = Kk.alpha   
                dtdt2 += ((Kk.alpha_prime_2(i, j) * alpha_k - Kk.alpha_prime(i) * Kk.alpha_prime(j)) / alpha_k ** 2 + Kk.beta_prime_2(i, j)) * self.dFda(k)
        
        return dtdt1 + dtdt2
    
def FillUpMatrices(F, ais, tis):
    F.Define(ais, tis)
    m = len(ais)
    n = len(tis)
    grad = np.array([F.dFda(i) for i in range(m)] + [F.dFdt(i) for i in range(m)])

    H = np.zeros((m + n, m + n))
    
    for i in range(m):
        for j in range(m):
            H[i, j] = F.d2Fda2(i, j)

    for i in range(m, m + n):
        for j in range(n):
            H[i, j] = F.d2Fdtda(i - m, j)
            H[j, i] = H[i, j]
            
    for i in range(m, m + n):
        for j in range(m, m + n):
            H[i, j] = F.d2Fdt2(i - m, j - m)
    #print(grad)
    #print("H", H)
    #print()
    #print("H_inv", np.linalg.inv(H))
    #print()
    #print("gradient: ", grad)
    #para
    return grad, np.linalg.inv(H)       
    
def GetExponentialsCoefficients(functional, a0, t0):
    tol = 1e-9
    max_iter = 100
    x = np.array(a0 + t0)
    x_old = np.array(a0 + t0) 
    
    still_changes = True
    iteration = 0

    while still_changes and iteration < max_iter:
        iteration += 1
        grad, H_inv = FillUpMatrices(functional, x[:len(a0)], x[len(a0):])
        x -= H_inv.dot(grad)    
        still_changes = np.linalg.norm(x - x_old) > tol
        x_old[:] = x[:]
        
    a0[:] = x[:len(a0)]
    t0[:] = x[len(a0):] 
    
# MAIN
#****************************************************************************************************************************************************************************************

def TimesToTaus(times):
    print(times)
    taus = [0.] * len(times)
    
    for i in range(len(times) - 1, 0, -1):
        taus[i] = math.log(times[i] - times[i - 1])

    print(times[0])
    taus[0] = math.log(times[0])
    return taus

def TausToTimes(taus):
    return [sum([math.exp(t) for t in taus[:i + 1]]) for i in range(len(taus))]

if __name__ == "__main__":
    tis = [0.1, 0.3, 1., 3., 5.,10., 40., 190., 1000., 6500., 50000.]
    a0 = [0.4 for t in tis]
    #a0 = [ 0.23477446,  0.28549392,  0.28479113, 0.3, 0.26149251,  0.32055117,  0.35351918, 0.3963018,   0.42237921,  0.48282255,  0.63471099]    
    tis = [0.1, 2.0, 20, 100]
    m = 5
    a0 = [1.0 for i in range(m)]  
    tis = [0.1 * math.exp(i) for i in range(m)]
    tis = TimesToTaus(tis)
    tol = 1e-9
    tol_residual = 1e-6
    max_iter = 30
    still_changing = True
    a = np.array(a0 + tis)
    a_old = np.array(a0 + tis) 
    a_best = np.array(a0 + tis)  
    iteration = 0
    F = Functional()
    F.Define(a[:len(a0)], a[len(a0):])
    print("calculating Fmod...")
    mod_residual = F.Fmod()
    print("calculating F...")
    best_residual = F.F()
    old_residual = best_residual
    print("RESIDUAL", best_residual)
    print("MOD_RESIDUAL", mod_residual)
    gamma_0 = 0.5
    
    while still_changing and iteration < max_iter:
        iteration += 1
        grad, H_inv = FillUpMatrices(F, a[:len(a0)], a[len(a0):])
        p = H_inv.dot(grad)   
        a -= p
        F.Define(a[:len(a0)], a[len(a0):])
        print("\ncalculating F...")
        residual =  F.F()
        print("calculating Fmod...\n")
        mod_residual = F.Fmod()
        
        if residual < best_residual:
            best_residual = residual
            a_best[:] = a[:]
        
        gradient_norm = sum([abs(float(g)) for g in grad])        
        print("\nGradient Norm", gradient_norm)
        print("RESIDUAL", residual)
        print("BEST_RESIDUAL", best_residual)
        print("MOD_RESIDUAL", mod_residual)
        print("\nBest times so far\n", TausToTimes(a_best[len(a0):]))
        print("\nBest coefficients so far\n", a_best[:len(a0)])
        print("\nCurrent times\n", TausToTimes(a[len(a0):]))
        print("\nCurrent coefficients\n", a[:len(a0)])
        still_changing = gradient_norm > tol 
        a_old[:] = a[:]
        
        
    #while still_changing and iteration < max_iter:
        #iteration += 1
        #grad, H_inv = FillUpMatrices(F, a[:len(a0)], a[len(a0):])
        #p = H_inv.dot(grad)          
        
        #not_improving = True
        #gamma = 1.
        #not_working = False
        #while not_improving:
            #a[:] = a_old[:] 
            #a -= gamma * p
            #F.Define(a[:len(a0)], a[len(a0):])
            #residual = F.F()            
            #not_improving = residual > old_residual
            #print(not_improving)
            #if residual < best_residual:
                #best_residual = residual
                #a_best[:] = a[:]
            #gamma *= 0.5
            #if gamma < 1e-8:
                #not_working = True
                #break
        #if not_working:
            #break
        #old_residual = residual        
        #mod_residual = F.Fmod()
        #print("\n RESIDUAL", residual)
        #print("BEST_RESIDUAL", best_residual)
        #print("MOD_RESIDUAL", mod_residual)
        ##improving = True
        ##gamma = - 1.0
        ##default_residual = residual
        ##while gamma < 2.0:   
            ##print("\ntrying to get closer...")
            ##print("gamma: ", gamma)
            ##print("current residual:", residual)
            ##print("best residual so far:", best_residual)
            ##a[:] = a_old[:] - gamma * p[:]
            ##gamma += gamma_0
            ##F.Define(a[:len(a0)], a[len(a0):])
            ##new_residual = F.F()
            ##if new_residual < 0:
                ##raise ValueError("You are getting negative residuals!!")
            ##improving = new_residual < residual or default_residual < residual
            ##if new_residual < best_residual:
                ##best_residual = residual
                ##a_best[:] = a[:]
            ##if best_residual < tol_residual:                    
                ##break
            ##residual = new_residual            

        #if best_residual < tol_residual:
            #print("The residual has become small enogh, quitting iterations")
            #break
        #still_changing = np.linalg.norm(a - a_old) > tol 
        #print("Change: ", np.linalg.norm(a - a_old))
        #print("best_residual_so_far", best_residual)
        #a_old[:] = a[:]

    print("\nBEST_RESIDUAL", best_residual)
    print("still changing: ", still_changing)    
    print("\nCOEFFICIENTS: ", a_best[:len(a0)])
    print(" \TIMES", TausToTimes(a_best[len(a0):]))
