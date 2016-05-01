import math 
import cmath
import mpmath
import matplotlib.pyplot as plt
from bigfloat import *
import numpy as np

def ExactIntegrationOfSinusKernel(t, t0 = 0):
    up_to_t  = 0.5 * math.pi * math.sqrt(t) * (mpmath.angerj(0.5, t) - mpmath.angerj(-0.5, t))
    up_to_t0 = 0.5 * math.pi * math.sqrt(t0) * (mpmath.angerj(0.5, t0) - mpmath.angerj(-0.5, t0))
    return up_to_t - up_to_t0

f = math.sin
final_time = 1.0
n_steps = 4
dt = final_time / n_steps
times = [dt * i for i in range(n_steps)]
times.append(final_time)
t_win = 0.3 * times[- 1]
for i in range(len(times)):
    if times[i] >= t_win:
        break
old_times = [time * times[i] / times[-1] for time in times]
recent_times = [times[i] + time * (times[-1] - times[0] - times[i]) / (times[-1] - times[0]) for time in times]
q = - 0.5
t = times[- 1]
a = times[0]
N = len(times) - 1
h = (t - a) / N
coeff = h ** (- q)
values = [gamma(k - q) / gamma(k + 1) * f(t - k * h) for k in range(N)]

t = recent_times[- 1]
a = recent_times[0]
N = len(recent_times) - 1
h = (t - a) / N
recent_coeff = h ** (- q)
recent_values = [gamma(k - q) / gamma(k + 1) * (f(t - k * h) - f(a)) for k in range(N)]
recent_correction = t ** (- q) / gamma(1 - q) * f(a)
#recent_correction = 2*sqrt(t-a)*f(a)#(t - a) ** (1 + q) / (1 + q) * f(a)

t = old_times[- 1]
a = old_times[0]
N = len(old_times) - 1
h = (t - a) / N
old_coeff = h ** (- q)
old_values = [gamma(k - q) / gamma(k + 1) * f(t - k * h) for k in range(N)]

print("times", times)
print("old_times", old_times)
print("recent_times", recent_times)
print("exact_old", ExactIntegrationOfSinusKernel(recent_times[0]))
print("exact_recent", ExactIntegrationOfSinusKernel(1, recent_times[0]))
print("exact all", ExactIntegrationOfSinusKernel(1.))
print("exact all minus recent", ExactIntegrationOfSinusKernel(1) - ExactIntegrationOfSinusKernel(recent_times[0]))
print("old", old_coeff * sum(old_values))
print("recent", recent_coeff * sum(recent_values) + recent_correction)
print("all", coeff * sum(values))
print("all minus old", coeff * sum(values) - old_coeff * sum(old_values))