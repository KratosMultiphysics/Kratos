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

best_as_L1 = []
best_ts_L1 = []
error_bounds_L1 = []
best_as_L1.append([1.055699152])
best_ts_L1.append([1.656571537])
error_bounds_L1.append(0.593414984)
best_as_L1.append([0.595936548, 0.765862627])
best_ts_L1.append([0.758737731, 8.130844515])
error_bounds_L1.append(0.274615576)
best_as_L1.append([0.457076294, 0.52049493, 0.730234918])
best_ts_L1.append([0.523735503, 3.463557465, 35.209010652])
error_bounds_L1.append(0.141643164)
best_as_L1.append([0.378123792, 0.420250984, 0.515234662, 0.747882647])
best_ts_L1.append([0.377168054, 1.883663548, 12.085534613, 127.522988007])
error_bounds_L1.append(0.077836293)
best_as_L1.append([0.377585589, 0.389837358, 0.414949491, 0.503856364, 0.607332741])
best_ts_L1.append([0.361079805, 1.758926107, 8.640539541, 57.10122954, 448.083463993])
error_bounds_L1.append(0.050354036)
best_as_L1.append([0.338300743, 0.345524197, 0.368960284, 0.368902685, 0.432603065, 0.771632072])
best_ts_L1.append([0.346126312, 1.386290002, 5.934710427, 27.706980453, 132.567423265, 1371.238854372])
error_bounds_L1.append(0.032416167)
best_as_L1.append([0.28596607, 0.296473585, 0.323588913, 0.360741831, 0.417056856, 0.514260513, 0.746597256])
best_ts_L1.append([0.222165932, 0.730394698, 2.617417995, 10.658764953, 52.200869695, 340.581473772, 3862.218173227])
error_bounds_L1.append(0.017434941)
best_as_L1.append([0.085504456, 0.242521351, 0.284508641, 0.501709401, 0.09440901, 0.570264153, 0.504139981, 0.851415185])
best_ts_L1.append([0.437648426, 0.451108658, 1.280099948, 6.252764672, 27.644724929, 140.107295013, 911.529174141, 10266.829839251])
error_bounds_L1.append(0.087813976)

hinsberg_as = [0.23477481312586, 0.28549576238194, 0.28479416718255, 0.26149775537574, 0.32056200511938, 0.35354490689146, 0.39635904496921, 0.42253908596514, 0.48317384225265, 0.63661146557001]
hinsberg_ts = [0.1, 0.3, 1., 3., 10., 40., 190., 1000., 6500., 50000.]


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
    print("m = ", m, "ERROR_L2: ", errors[i])

numbers_L1 = []
errors_L1 = []
window_kernel_errors = [exact_integral for a in best_as_L1]
for i in range(len(best_as_L1)):
    ais = best_as_L1[i]
    tis = best_ts_L1[i]
    tis = [ti * t_win for ti in tis]
    m = len(ais)
    approximate_inegral = float(IntegrationWithSumOfExponentialsKernel(t, t0, t_tail, ais, tis))
    errors_L1.append(abs(exact_integral - approximate_inegral) * math.sqrt(t_win))
    numbers_L1.append(m)
    print("m = ", m, "ERROR_L1: ", errors_L1[i])

approximate_inegral_hinsberg = float(IntegrationWithSumOfExponentialsKernel(t, t0, t_tail, hinsberg_as, hinsberg_ts))
error_hinsberg = abs(exact_integral - approximate_inegral_hinsberg) * math.sqrt(t_win)

plt.plot(numbers, errors, label = "error")
plt.plot(numbers, error_bounds, label = "bound")
plt.plot(numbers, window_kernel_errors, label = "window kernel error")
plt.plot(numbers_L1, errors_L1, label = "error_L1")
plt.plot(numbers_L1, error_bounds_L1, label = "bound_L1")
plt.plot([10], [error_hinsberg], 'bo',  label = "Hinsberg")
plt.xlim([0, 11])
plt.semilogy()
plt.legend()
plt.show()