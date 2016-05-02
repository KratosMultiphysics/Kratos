import math 
import cmath
import mpmath
import matplotlib.pyplot as plt
from bigfloat import *
import numpy as np

def ExactIntegrationOfSinusKernel(t, a = None, b = None):
    with precision(300):
        if a == None and b == None:
            return 0.5 * math.pi * math.sqrt(t) * (mpmath.angerj(0.5, t) - mpmath.angerj(-0.5, t))
        if a == None and b != None:
            a = t
        elif b == None:
            b = 0.
        mpmath.mp.pretty = True
        pi = mpmath.mp.pi
        pi = +pi
        fcos = mpmath.fresnelc
        fsin = mpmath.fresnels
        
        arg_a = mpmath.sqrt(2 * (t - a) / mpmath.mp.pi)
        arg_b = mpmath.sqrt(2 * (t - b) / mpmath.mp.pi) 
        return mpmath.sqrt(2 * mpmath.mp.pi) * ((fsin(arg_b) - fsin(arg_a)) * mpmath.cos(t) + (fcos(arg_a) - fcos(arg_b)) * mpmath.sin(t))
f = math.sin
final_time = 1.0
n_steps = 100
dt = final_time / n_steps
times = [dt * i for i in range(n_steps)]
times.append(final_time)
t_win = times[0] + 0.5 * times[- 1]
for i in range(len(times)):
    if times[i] >= t_win:
        break
old_times = [time * times[i] / times[-1] for time in times]
recent_times = [times[i] + time * (times[-1] - times[0] - times[i]) / (times[-1] - times[0]) for time in times]
q = - 0.5
t = times[- 1]
N = len(times) - 1
h = t / N
coeff = h ** (- q)
values = [gamma(k - q) / gamma(k + 1) * f(t - k * h) for k in range(N)]

t = old_times[- 1]
h = t / N
old_coeff = h ** (- q)
old_values = [gamma(k - q) / gamma(k + 1) * (f(t - k * h) - f(t)) for k in range(N + 1)]
old_correction = 2 * (sqrt(final_time) - sqrt(final_time - t)) * f(t)

t = recent_times[- 1]
a = recent_times[0]
t -= a
h = t / N
initial_approx_deriv = cos(a)
constant_initial_correction = 2 * sqrt(t) * f(t)#t ** (- q) / gamma(1 - q) * f(a)
linear_initial_correction = 2. / 3 * sqrt(t) * (2 * t) * initial_approx_deriv 
#linear_initial_correction = t ** (1 - q) / gamma(2 - q) * initial_approx_deriv
linear_correction_option = 0

recent_coeff_1 = h ** (- q)
recent_values_1 = [gamma(k - q) / gamma(k + 1) * (f(t - k * h + a) - f(t)) for k in range(N + 1)]
#recent_correction = gamma(- q) * t ** (- q) / gamma(1 - q) * f(a)
recent_correction_1 = constant_initial_correction

recent_coeff_2 = h ** (- q) * gamma(- q)
recent_values_2 = [(- 1) ** k * gamma(q + 1) / (gamma(k + 1) * gamma(q - k + 1)) * (f(t - (k  - 0.5 * q) * h + a) - f(a) - linear_correction_option * (t - (k  - 0.5 * q) * h) * initial_approx_deriv) for k in range(N)]
#recent_correction = gamma(- q) * t ** (- q) / gamma(1 - q) * f(a)
recent_correction_2 = initial_value_correction =  gamma(- q) * (constant_initial_correction + linear_correction_option * linear_initial_correction)

final_time = times[- 1]
t_win = recent_times[0]

print("times", times)
print("old_times", old_times)
print("recent_times", recent_times)
print("EXACT_old", ExactIntegrationOfSinusKernel(final_time, 0., t_win))
print("EXACT_recent", ExactIntegrationOfSinusKernel(final_time, t_win, final_time))
print("EXACT all", ExactIntegrationOfSinusKernel(final_time))
print("old", old_coeff * sum(old_values) + old_correction)
print("RECENT_1", recent_coeff_1 * sum(recent_values_1) + recent_correction_1)
print("RECENT_2", recent_coeff_2 * sum(recent_values_2) + recent_correction_2)
print("all", coeff * sum(values))
print("all minus old", coeff * sum(values) - old_coeff * sum(old_values))