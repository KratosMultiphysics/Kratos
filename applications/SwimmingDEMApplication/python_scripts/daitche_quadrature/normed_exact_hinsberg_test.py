import math
import cmath
import mpmath
import matplotlib.pyplot as plt
import quadrature as quad
import optimal_points as op
import scipy as sp
from scipy.special import lambertw
from scipy.optimize import fsolve
from matplotlib import rc

class ProblemParameters:
    def __init__(self):
        pass


def bisection(f, a, b, tol):
    assert f(a) * f(b) < 0.0
    c = 0.5 * (a + b)
    it = 0
    while 0.5 * (b - a) > tol and it < 100:
        if f(c) == 0:
            return c
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c
        c = 0.5 * (a + b)
        it += 1

    return c

def CalculateErrors(points_set, pp):
    ExactIntegral = quad.ExactIntegrationOfSinus(pp.end_time, pp.initial_time)
    pp.lower_limits_for_tail = [0.0] * len(points_set.AsTs)

    for m in points_set.exponential_indices:
        tis = [points_set.AsTs[m][1][i] for i in range(len(points_set.AsTs[m][0]))]
        t_min = tis[0]
        t_max = tis[- 1]

        def K_prime(t):
            ais = [points_set.AsTs[m][0][i] for i in range(len(points_set.AsTs[m][0]))]
            alphas = [math.sqrt(math.exp(1) / ti) for ti in tis]
            betas = [- 0.5 / ti for ti in tis]
            return sum([ais[i] * alphas[i] * betas[i] * math.exp(betas[i] * t) for i in range(len(points_set.AsTs[m][0]))])
        rho = 2.
        t = t_max

        while rho > 1:
            t += t
            value_at_last_tg_point_K_prime = K_prime(t)
            value_at_last_tg_point_KB_prime = - 1.0 / 2 / t ** 1.5
            rho = abs(value_at_last_tg_point_K_prime / value_at_last_tg_point_KB_prime)

        security_coeff = 2.0

        pp.lower_limits_for_tail[m] = security_coeff * t
        #print('rho_calc', math.exp(beta_last_exp*(pp.lower_limits_for_tail[m] - t_max)))
        #print('pp.lower_limits_for_tail[m]', pp.lower_limits_for_tail[m])
        #print('lets see', 'm = ', m + 1, 'K / K_B (t) = ', abs(K_prime(pp.lower_limits_for_tail[m]) * 2 * pp.lower_limits_for_tail[m] ** 1.5))

    for m in points_set.exponential_indices:
        F_w = quad.ExactIntegrationOfSinus(pp.end_time, pp.end_time_minus_tw)
        tis = [pp.t_w * ti for ti in points_set.AsTs[m][1]]
        F_tail = quad.ExactIntegrationOfTail(pp.end_time, pp.end_time_minus_tw, pp.initial_time, points_set.AsTs[m][0], tis)
        ApproxIntegral = float(F_w) + float(F_tail)
        error_bound = ObjectiveFunction(points_set, m, pp)
        Error = abs(ExactIntegral - ApproxIntegral)
        points_set.Errors.append(Error)
        points_set.ErrorBounds.append(error_bound)

def ExactErrorTail(points_set, m_index, lower_limit):
    ais = [points_set.AsTs[m_index][0][i] for i in range(len(points_set.AsTs[m_index][0]))]
    tis = [points_set.AsTs[m_index][1][i] for i in range(len(points_set.AsTs[m_index][0]))]
    alphas = [math.sqrt(math.exp(1) / ti) for ti in tis]
    betas = [- 0.5 / ti for ti in tis]
    return 1. / math.sqrt(lower_limit) - sum([ais[i] * alphas[i] * math.exp(betas[i] * lower_limit) for i in range(len(points_set.AsTs[m_index][0]))])

def K(points_set, m_index, t, pp):
        ais = [points_set.AsTs[m_index][0][i] for i in range(len(points_set.AsTs[m_index][0]))]
        tis = [points_set.AsTs[m_index][1][i] for i in range(len(points_set.AsTs[m_index][0]))]
        alphas = [math.sqrt(math.exp(1) / ti) for ti in tis]
        betas = [- 0.5 / ti for ti in tis]
        return sum([ais[i] * alphas[i] * math.exp(betas[i] * t) for i in range(len(points_set.AsTs[m_index][0]))])

def GetMaxTime(max_a, min_b, m, pp):
        from scipy.special import lambertw
        t_intersection = (2 * min_b * pp.end_time - lambertw(2 * min_b / ((m + 1) * min_a) ** 2)) / (2 * min_b)
        return t_intersection

def PerformQuadratureOfObjectiveFunctionSecondTerm(points_set, m_index, a, b):

    def DObjectiveFunctionDt(t):
        ais = [points_set.AsTs[m_index][0][i] for i in range(len(points_set.AsTs[m_index][0]))]
        tis = [points_set.AsTs[m_index][1][i] for i in range(len(points_set.AsTs[m_index][0]))]
        alphas = [math.sqrt(math.exp(1) / ti) for ti in tis]
        betas = [- 0.5 / ti for ti in tis]
        return abs(- 0.5 * t ** (- 1.5) - sum([ais[i] * alphas[i] * betas[i] * math.exp(betas[i] * t) for i in range(len(points_set.AsTs[m_index][0]))]))

    from scipy.integrate import quad
    value, error = quad(DObjectiveFunctionDt, a, b, limit = 500000)
    return value, error

def ObjectiveFunction(points_set, m_index, pp):
    return 1.0 / math.sqrt(pp.t_w) * (abs(1 - K(points_set, m_index, 1, pp)) + PerformQuadratureOfObjectiveFunctionSecondTerm(points_set, m_index, 1, pp.lower_limits_for_tail[m_index])[0] + ExactErrorTail(points_set, m_index, pp.lower_limits_for_tail[m_index]))

def Average(vector_container, current_values, k_calc, pp):
    i_tw = k_calc - 1
    current_vector_length = len(current_values)

    if len(vector_container[i_tw]) != current_vector_length:
        vector_container[i_tw] = [0] * current_vector_length

    for i in range(current_vector_length):
        old_value = float(vector_container[i_tw][i])
        if pp.error_norm_type == 'L1':
            vector_container[i_tw][i] = 1.0 / k_calc * ((i_tw) * vector_container[i_tw][i] + current_values[i])
        elif pp.error_norm_type == 'max':
            vector_container[i_tw][i] = max(abs(vector_container[i_tw][i]), abs(current_values[i]))
        pp.rate_of_change += abs((vector_container[i_tw][i] - old_value) / (0.5 * (old_value + vector_container[i_tw][i])))

class HinsbergPointsSetGivenNorm:
    def __init__(self, code):
        self.code = code
        self.AsTs = []
        self.Errors = []
        self.ErrorBounds = []

pp = ProblemParameters()
pp.n_samples = 40
pp.t_w_min = 0.00001
pp.n_doublings = 6
pp.error_norm_type = 'max'
pp.initial_number_of_periods = 1
phases = [i / pp.n_samples * 2 * math.pi  for i in range(pp.n_samples)]
t_ws = [pp.t_w_min * 10 ** k for k in range(pp.n_doublings)]
#t_ws = [pp.t_w_min]
exponential_numbers_t_norm = [0]
exponential_numbers_abs_norm = [0]
exponential_numbers_hinsberg_norm = [0]
norm_of_errors_t_norm = [[]] * len(t_ws)
norm_of_errors_abs_norm = [[]] * len(t_ws)
norm_of_errors_hinsberg_norm = [[]] * len(t_ws)
norm_of_bounds_t_norm = [[]] * len(t_ws)
norm_of_bounds_abs_norm = [[]] * len(t_ws)
norm_of_bounds_hinsberg_norm = [[]] * len(t_ws)
pp.rate_of_change = 0.0
k_sample = 0

for phase in phases:
    print()
    print('COMPUTING SAMPLE ' + str(k_sample + 1) + ' OF ' + str(len(phases)) + '...')
    print('----------------------------------------------')
    pp.end_time = pp.initial_number_of_periods * 2 * math.pi + phase
    pp.initial_time = 'MinusInfinity'
    k_calc = 1
    k_max = len(t_ws)

    for t_w in t_ws:
        print('Calculation ' + str(k_calc) + ' of ' + str(k_max) + '.')
        print('t_w = ', t_w)
        pp.t_w = t_w
        pp.end_time_minus_tw = pp.end_time - pp.t_w

        t_norm_set = HinsbergPointsSetGivenNorm('t_norm')
        abs_norm_set = HinsbergPointsSetGivenNorm('abs_norm')
        hinsberg_set = HinsbergPointsSetGivenNorm('hinsberg_norm')
        op.FillPoints(t_norm_set, t_norm_set.code)
        op.FillPoints(abs_norm_set, abs_norm_set.code)
        op.FillPoints(hinsberg_set, hinsberg_set.code)
        t_norm_set.exponential_indices = [m for m in range(len(t_norm_set.AsTs))]
        abs_norm_set.exponential_indices = [m for m in range(len(abs_norm_set.AsTs))]
        hinsberg_set.exponential_indices = [0]

        CalculateErrors(t_norm_set, pp)
        CalculateErrors(abs_norm_set, pp)
        CalculateErrors(hinsberg_set, pp)
        t_norm_set.exponential_numbers = [m + 1 for m in t_norm_set.exponential_indices]
        abs_norm_set.exponential_numbers = [m + 1 for m in abs_norm_set.exponential_indices]
        hinsberg_set.exponential_numbers = [10]
        exponential_numbers_t_norm        = t_norm_set.exponential_numbers
        exponential_numbers_abs_norm      = abs_norm_set.exponential_numbers
        exponential_numbers_hinsberg_norm = hinsberg_set.exponential_numbers
        Average(norm_of_errors_abs_norm, abs_norm_set.Errors, k_calc, pp)
        Average(norm_of_errors_t_norm, t_norm_set.Errors, k_calc, pp)
        Average(norm_of_errors_hinsberg_norm, hinsberg_set.Errors, k_calc, pp)
        Average(norm_of_bounds_abs_norm, abs_norm_set.ErrorBounds, k_calc, pp)
        Average(norm_of_bounds_t_norm, t_norm_set.ErrorBounds, k_calc, pp)
        Average(norm_of_bounds_hinsberg_norm, hinsberg_set.ErrorBounds, k_calc, pp)
        k_calc += 1
    k_sample += 1
    print('----------------------------------------------')
    print('Rate of change: ' + str(pp.rate_of_change))
    pp.rate_of_change = 0.0

k_calc = 1
for t_w in t_ws:
    size_factor = k_max
    line_width = (size_factor * k_calc / k_max) ** 1.2
    small_marker_size = 4 + 2 * (size_factor * k_calc / k_max) ** 1.2
    big_marker_size = 4 + 3 * (size_factor * k_calc / k_max) ** 1.2
    maker_width = 0.75 * (size_factor * k_calc / k_max) ** 1.2
    if k_calc == 3:
        plt.plot(exponential_numbers_t_norm, norm_of_errors_t_norm[k_calc - 1], markersize = small_marker_size, linewidth = line_width, color='b', marker = 'o', label= '$I_{2t}$')
        plt.plot(exponential_numbers_abs_norm, norm_of_errors_abs_norm[k_calc - 1], markersize = big_marker_size, linewidth = line_width, color='k', marker = '*', label= '$I_1$')
        plt.plot(exponential_numbers_hinsberg_norm, norm_of_errors_hinsberg_norm[k_calc - 1], markersize = 20, mew = maker_width, color='g', marker = 'x', label= '$I_{2tH}$, points by van Hinsberg et al. ')

        plt.plot(exponential_numbers_t_norm, norm_of_bounds_t_norm[k_calc - 1], color='b', linestyle='--', label= r'$I_{2t}$ (bound)')
        plt.plot(exponential_numbers_abs_norm, norm_of_bounds_abs_norm[k_calc - 1], color='k', linestyle=':', label= '$I_1$ (bound)')
        plt.plot(exponential_numbers_hinsberg_norm, norm_of_bounds_hinsberg_norm[k_calc - 1], markersize = 10, color='g', linestyle='--', marker = 'v', label= '$I_{2tH}$, points by van Hinsberg et al. (bound)')
    else:
        plt.plot(exponential_numbers_t_norm, norm_of_errors_t_norm[k_calc - 1], markersize = small_marker_size, linewidth = line_width, color='b', marker = 'o')
        plt.plot(exponential_numbers_abs_norm, norm_of_errors_abs_norm[k_calc - 1], markersize = big_marker_size, linewidth = line_width, color='k', marker = '*')
        plt.plot(exponential_numbers_hinsberg_norm, norm_of_errors_hinsberg_norm[k_calc - 1], markersize = 20, mew = maker_width, color='g', marker = 'x')

        plt.plot(exponential_numbers_t_norm, norm_of_bounds_t_norm[k_calc - 1], color='b', linestyle='--')
        plt.plot(exponential_numbers_abs_norm, norm_of_bounds_abs_norm[k_calc - 1], color='k', linestyle=':')
        plt.plot(exponential_numbers_hinsberg_norm, norm_of_bounds_hinsberg_norm[k_calc - 1], markersize = 10, color='g', linestyle='--', marker = 'v')
    k_calc += 1

plt.xlim(0, 11)
plt.ylabel('$E$', fontsize=40, labelpad=25)
plt.xlabel('$m$', fontsize=40, labelpad=5)
plt.tick_params(axis='both', which='major', labelsize=35)
plt.semilogy()
plt.legend(loc = 'lower left', prop={'size':22.5},frameon=False)
ax = plt.gca()
ax.tick_params(axis='x', pad=20)
ax.tick_params(axis='y', pad=10)
plt.savefig('exact_error_with_sinus_min_t_w=' + str(pp.t_w_min) + '_max_tw=' + str(pp.t_w) + '.eps', format='eps', dpi=1000)
plt.show()
