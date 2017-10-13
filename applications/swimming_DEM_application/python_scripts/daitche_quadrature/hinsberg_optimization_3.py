import math 
import cmath
import mpmath
import matplotlib.pyplot as plt
from bigfloat import *
import numpy as np
import scipy
import cmath

def DawsonIntegral(x):
    return 0.5 * sqrt(math.pi) * exp(- x ** 2) * scipy.special.erfi(x)

class K_B:
    def __init__(self):
        pass
    def f(self, t):
        return t ** (- 0.5)
    
    def df(self, t):
        return - 0.5 * t ** (- 1.5)

class K_component:
    def __init__(self, taui, ai):
        self.alpha = exp(0.5 * (1 - taui))
        self.alpha_prime = - 0.5 * self.alpha
        self.alpha_prime_2 = 0.25 * self.alpha
        self.beta = - 0.5 * exp(- taui)
        self.beta_prime = - self.beta
        self.beta_prime_2 = self.beta
        self.a = ai
        
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
        self.Ks = [K_component(tis[i], ais[i]) for i in range(len(tis))]
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
        Ki = K.Ks[i]
        alpha_i = Ki.alpha   
        alpha_prime_i = Ki.alpha_prime
        beta_prime_i = Ki.beta_prime
        a_i = Ki.a
        return a_i * (alpha_prime_i / alpha_i + beta_prime_i) * self.dFda(i)
    
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
        Ki = K.Ks[i]   
        alpha_i = Ki.alpha           
        alpha_prime_i = Ki.alpha_prime
        beta_prime_i = Ki.beta_prime        
        a_i = Ki.a        
        dtda = a_i *  self.d2Fda2(i, j)
        if i == j:
            dtda +=  self.dFda(i)
            
        dtda *= (alpha_prime_i / alpha_i + beta_prime_i)
        return dtda

    def d2Fdt2(self, i, j):
        K = self.K
        Ki = K.Ks[i]      
        alpha_i = Ki.alpha           
        alpha_prime_i = Ki.alpha_prime
        beta_prime_i = Ki.beta_prime
        alpha_prime_2_i = Ki.alpha_prime_2
        beta_prime_2_i = Ki.beta_prime_2 
        a_i = Ki.a               
        dtdt =  self.d2Fdtda(j, i)
        dtdt *= (alpha_prime_i / alpha_i + beta_prime_i) * a_i
        if i == j:
            dtdt += a_i * ((alpha_prime_2_i * alpha_i - alpha_prime_i ** 2) / alpha_i ** 2 + beta_prime_2_i) *  self.dFda(i)
        return dtdt

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
      
    print("H", H)
    print()
    print("H_inv", np.linalg.inv(H))
    print()
    print("gradient: ", grad)
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
if __name__ == "__main__":
    tis = [0.1, 0.3, 1., 3., 10., 40., 190., 1000., 6500., 50000.]
    a0 = [ 0.23477446,  0.28549392,  0.28479113,  0.26149251,  0.32055117,  0.35351918, 0.3963018,   0.42237921,  0.48282255,  0.63471099]    
    tis = [0.01, 0.02, 0.03]
    tis = [math.log(t) for t in tis]
    a0 = [0.1,1.0,10.0]
    tol = 1e-9
    max_iter = 50
    still_changes = True
    a = np.array(a0 + tis)
    a_old = np.array(a0 + tis) 
    iteration = 0
    F = Functional()

    while still_changes and iteration < max_iter:
        iteration += 1
        grad, H_inv = FillUpMatrices(F, a[:len(a0)], a[len(a0):])
        a -= H_inv.dot(grad)  
        still_changes = np.linalg.norm(a - a_old) > tol 
        print("Change: ", np.linalg.norm(a - a_old))
        a_old[:] = a[:]


    print("a coefficients: ", a[:len(a0)])
    print("times: ", [math.exp(tau) for tau in a[len(a0):]])
    print("still changing: ", still_changes)