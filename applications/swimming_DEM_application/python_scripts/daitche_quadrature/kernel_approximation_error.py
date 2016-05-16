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

def ExactIntegrationOfSinus(t, a = None, b = None):
    with precision(300):
        if a == None and b == None:
            return 0.5 * math.pi * math.sqrt(t) * (mpmath.angerj(0.5, t) - mpmath.angerj(- 0.5, t))
        elif a == None and b != None:
            a = 0
        elif b == None:
            b = t
        mpmath.mp.dps = 50
        mpmath.mp.pretty = True
        pi = mpmath.mp.pi
        pi = +pi
        fcos = mpmath.fresnelc
        fsin = mpmath.fresnels        
        arg_b = mpmath.sqrt(2 * (t - b) / pi) 
        
        if a == "MinusInfinity":
            return mpmath.sqrt(0.5 * mpmath.mp.pi) * (- 2 * mpmath.sin(t) * fcos(arg_b) + 2 * mpmath.cos(t) * fsin(arg_b) + mpmath.sin(t) - mpmath.cos(t))
        else:
            arg_a = mpmath.sqrt(2 * (t - a) / pi)
            return mpmath.sqrt(2 * mpmath.mp.pi) * ((fsin(arg_b) - fsin(arg_a)) * mpmath.cos(t) + (fcos(arg_a) - fcos(arg_b)) * mpmath.sin(t))      

def ExactIntegrationOfSinusWithExponentialKernel(ti, t, a, b):
    alpha = sqrt(exp(1) / ti)
    beta = - 0.5 / ti
    exp_b = exp(beta * (t - b))
    sb = sin(b)
    cb = cos(b)
    coeff = 1. / (beta ** 2 + 1)    
    
    if a == "MinusInfinity":
        return - alpha * coeff * exp_b * (beta * sb + cb)
    else:
        exp_a = exp(beta * (t - a))
        sa = sin(a)
        ca = cos(a)
        return alpha * coeff * (exp_a * (ca - beta * sa) - exp_b * (beta * sb - cb))

def IntegrationWithSumOfExponentialsKernel(t, a, b, ais, tis):
    integral = 0.0
    for i in range(len(tis)):
        ti = tis[i]
        ai = ais[i]
        integral += ai * ExactIntegrationOfSinusWithExponentialKernel(ti, t, a, b)
        
    return integral

t = 5.0
t_win = 3.0
t0 = "MinusInfinity"
t_tail = t - t_win

best_as = []
best_ts = []
error_bounds = []
best_as.append([0.93847245053])
best_ts.append([1.43003412613])
error_bounds.append(0.630184409394704)
best_as.append([0.54705976442, 0.84497676080])
best_ts.append([0.66668356676, 8.34248846822])
error_bounds.append(0.284141395048289)
best_as.append([0.43079707202, 0.53194030882, 0.80464753089])
best_ts.append([0.45214632116, 3.05971221535, 36.76947294003])
error_bounds.append(0.145887389065311)
best_as.append([0.37140492497, 0.42213046005, 0.52488268154, 0.78143202593])
best_ts.append([0.35050516570, 1.75257036603, 11.65284347454, 136.88606867688])
error_bounds.append(0.088736294760644)
best_as.append([0.3335748221, 0.3629345041, 0.4197267406, 0.5202038586, 0.7661065349])
best_ts.append([0.2904630745, 1.2037066711, 5.9371428024, 39.1460800428, 452.8400756807])
error_bounds.append(0.062538312240541)
best_as.append([0.3065996056, 0.3243562296, 0.3616042093, 0.4181372918, 0.5168269107, 0.7551323048])
best_ts.append([0.2504415371, 0.9103723648, 3.7209169239, 18.2757857451, 119.7884523811, 1370.3543385994])
error_bounds.append(0.051954819317917)
best_as.append([0.286024225, 0.2965345272, 0.3232766788, 0.360869665, 0.4169241996, 0.5142197299, 0.7467744804])
best_ts.append([0.2216550529, 0.7303660627, 2.6151175192, 10.6587706724, 52.2008522707, 340.5814361913, 3862.2182275084])
error_bounds.append(0.045539886788637)
best_as.append([0.2696042712, 0.275174292, 0.2954251215, 0.3228225133, 0.3602769943, 0.4159412411, 0.5121587415, 0.7402444628])
best_ts.append([0.1998261277, 0.6095050731, 1.9749864724, 7.0542458787, 28.7011577713, 140.2284577358, 911.5441865997, 10266.8313490027])
error_bounds.append(0.042541670257497 )

numbers = []
errors = []
exact_integral = float(ExactIntegrationOfSinus(t, t0, t_tail))

for i in range(len(best_as)):
    ais = best_as[i]
    tis = best_ts[i]
    tis = [ti * t_win for ti in tis]
    m = len(ais)
    approximate_inegral = float(IntegrationWithSumOfExponentialsKernel(t, t0, t_tail, ais, tis))
    errors.append(abs(exact_integral - approximate_inegral) * math.sqrt(t_win))
    numbers.append(m)
    print("m = ", m, "ERROR: ", errors[i])

plt.plot(numbers, errors, label = "error")
plt.plot(numbers, error_bounds, label = "bound")
plt.semilogy()
plt.legend()
plt.show()