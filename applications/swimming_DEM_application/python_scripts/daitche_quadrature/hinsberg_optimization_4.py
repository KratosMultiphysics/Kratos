import math 
import cmath
import mpmath
import matplotlib.pyplot as plt
from bigfloat import *
import numpy as np
import scipy
import scipy.optimize as opt
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

class Functional:
    def __init__(self):
        self.sqrt_pi = sqrt(math.pi)
    
    def Define(self, ais, tis):
        self.K = K_sum(ais, tis)
        self.K_1 = self.K.f(1)
        
    def F(self):
        K = self.K
        def f(t):
            return float(- 1. / t ** 1.5 - sum([K_comp.a * K_comp.df(t) for K_comp in K.Ks]))

        try:
            t_cross = FindZero(f, 2.0) 
            first_bit = - 2 + 2 / sqrt(t_cross) - sum([K_comp.a * (K_comp.df(t_cross) - K_comp.df(1.)) / K_comp.beta for K_comp in K.Ks])
            second_bit = 2 / sqrt(t_cross) + sum([- K_comp.a * K_comp.df(t_cross) / K_comp.beta for K_comp in K.Ks])            
            print("\nT_CROSS",t_cross)
            print("BEFORE",f(t_cross-0.5*t_cross))
            print("AFT",f(t_cross+0.5*t_cross))
            print("first bit",first_bit)
            print("second bit", second_bit)
            return float(first_bit + second_bit)
        except:
            return float(2 + sum([- K_comp.a * K_comp.df(1.) / K_comp.beta for K_comp in K.Ks]))
        
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
            dfdt += a_k * (Kk.alpha_prime(i) / alpha_k + Kk.beta_prime(i))
        return dfdt * self.dFda(i)
    
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
        
        for k in range(i + 1):
            Kk = K.Ks[k]
            alpha_k = Kk.alpha   
            a_k = Kk.a
            dtda1 += a_k * (Kk.alpha_prime(i) / alpha_k + Kk.beta_prime(i))
        dtda1 *= self.d2Fda2(i, j)
        dtda2 = 0.
        
        if i >= j:
            Kj = K.Ks[j]
            alpha_j = Kj.alpha   
            dtda2 += (Kj.alpha_prime(i) / alpha_j + Kj.beta_prime(i))
        
        dtda2 *= self.dFda(i)
        return dtda1 + dtda2

    def d2Fdt2(self, i, j):
        K = self.K
        dtdt1 = 0.
        
        for k in range(i + 1):
            Kk = K.Ks[k]
            alpha_k = Kk.alpha   
            a_k = Kk.a
            dtdt1 += a_k * (Kk.alpha_prime(i) / alpha_k + Kk.beta_prime(i))     
        dtdt1 *= self.d2Fdtda(j, i)
        dtdt2 = 0.
        if i >= j:
            for k in range(i + 1):
                Kk = K.Ks[k]
                alpha_k = Kk.alpha   
                dtdt2 += ((Kk.alpha_prime_2(i, j) * alpha_k - Kk.alpha_prime(i) * Kk.alpha_prime(j)) / alpha_k ** 2 + Kk.beta_prime_2(i, j))
            dtdt2 * self.dFdt(i)
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
    #tis = [0.1, 0.3, 1., 3., 10., 40., 190., 1000., 6500., 50000.]
    #a0 = [ 0.23477446,  0.28549392,  0.28479113,  0.26149251,  0.32055117,  0.35351918, 0.3963018,   0.42237921,  0.48282255,  0.63471099]    
    tis = [0.1, 1.0, 100]
    a0 = [1.,1.,1.]    
    tis = TimesToTaus(tis)
    tol = 1e-9
    max_iter = 100
    still_changes = True
    a = np.array(a0 + tis)
    a_old = np.array(a0 + tis) 
    iteration = 0
    F = Functional()

    while still_changes and iteration < max_iter:
        iteration += 1
        grad, H_inv = FillUpMatrices(F, a[:len(a0)], a[len(a0):])
        p = H_inv.dot(grad)  
        gamma = 0.1
        
        residual = F.F()
        print("RESIDUAL", residual)
        improving = True
        a -= p
        #while improving:
            #a[:] = a_old[:] - gamma * p[:]
            #F.Define(a[:len(a0)], a[len(a0):])
            #new_residual = F.F()
            #improving = abs(new_residual) < abs(residual)
            #residual = new_residual
            #print("RESIDUAL", residual)
            #gamma *= (1 + 0.5)
            
        still_changes = np.linalg.norm(a - a_old) > tol 
        print("Change: ", np.linalg.norm(a - a_old))
        a_old[:] = a[:]


    print("a coefficients: ", a[:len(a0)])
    print("times: ", TausToTimes(a[:len(a0)]))
    print("still changing: ", still_changes)