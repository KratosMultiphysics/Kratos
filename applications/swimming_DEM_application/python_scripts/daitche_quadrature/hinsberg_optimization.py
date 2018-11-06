import math 
import cmath
import mpmath
import matplotlib.pyplot as plt
from bigfloat import *
import numpy as np

class K_B:
    def __init__(self):
        pass
    def f(self, t):
        return t ** (- 0.5)
    
    def df(self, t):
        return - 0.5 * t ** (- 1.5)

class K_component:
    def __init__(self, ti, ai):
        self.alpha = sqrt(exp(1) / ti)
        self.beta = - 1. / (2 * ti)
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
    def __init__(self, tis, ais):
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
    def __init__(self, tis):
        self.tis = tis
        self.sqrt_pi = sqrt(math.pi)
    
    def Define(self, ais):
        self.K = K_sum(self.tis, ais)
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
        K_1 = self.K_1
        sqrt_pi = self.sqrt_pi        
        Ki = K.Ks[i]
        Ki_1 = Ki.f(1)
        alpha_i = Ki.alpha        
        beta_i = Ki.beta        
        sqrt_minus_betai = sqrt(- beta_i)
        t_Ki_prime_K_prime = 2 * beta_i * Ki_1 * sum([Kk.a * Kk.f(1.) * Kk.beta * (1. - Kk.beta - beta_i) / (Kk.beta + beta_i) ** 2 for Kk in K.Ks])
        
        return - 2 * (1 - K_1) * Ki_1 + alpha_i * beta_i * sqrt_pi * erfc(sqrt_minus_betai) / sqrt_minus_betai + t_Ki_prime_K_prime
    
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

def FillUpMatrices(F, a):
    F.Define(a)
    m = len(a)
    grad = np.array([F.dFda(i) for i in range(m)])
    
    H = np.zeros((m, m))
    
    for i in range(m):
        for j in range(m):
            H[i,j] = F.d2Fda2(i, j)

    return grad, np.linalg.inv(H)       
    
def GetExponentialsCoefficients(functional, a0):
    tol = 1e-12
    max_iter = 10
    a = np.array(a0)
    a_old = np.array(a0) 
    
    still_changes = True
    iteration = 0

    while still_changes and iteration < max_iter:
        iteration += 1
        grad, H_inv = FillUpMatrices(functional, a)
        a -= H_inv.dot(grad)    
        still_changes = np.linalg.norm(a - a_old) > tol
        a_old[:] = a[:]
        
    a0[:] = a[:]
    
# MAIN
#****************************************************************************************************************************************************************************************
if __name__ == "__main__":
    tis = [0.1, 0.3, 1., 3., 10., 40., 190., 1000., 6500., 50000.]
    a0 = [0.2 for ti in tis]
    tol = 1e-9
    max_iter = 10
    still_changes = True
    a = np.array(a0)
    a_old = np.array(a0) 
    iteration = 0
    F = Functional(tis)

    while still_changes and iteration < max_iter:
        iteration += 1
        grad, H_inv = FillUpMatrices(F, a)
        a -= H_inv.dot(grad)    
        still_changes = np.linalg.norm(a - a_old) > tol
        a_old[:] = a[:]
    print("a coefficients: ", a)
    print("still changing: ", still_changes)