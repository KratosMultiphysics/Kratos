import math
import cmath
import mpmath
import matplotlib.pyplot as plt
import quadrature as quad
import optimal_points as op
import scipy as sp

class ProblemParameters:
    def __init__(self):
        pass
    
def CalculateErrors(points_set, pp):
    ExactIntegral = quad.ExactIntegrationOfSinus(pp.final_time)

    for m in points_set.exponential_indices:
        F_w = quad.ExactIntegrationOfSinus(pp.final_time, pp.final_time_minus_tw)
        tis = [pp.t_w * ti for ti in points_set.AsTs[m][1]]
        F_tail = quad.ExactIntegrationOfTail(pp.final_time, pp.final_time_minus_tw, 0, points_set.AsTs[m][0], tis)
        ApproxIntegral = float(F_w) + float(F_tail)
        error_bound = ObjectiveFunction(points_set, m, pp)
        Error = abs((ExactIntegral - ApproxIntegral))
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

def PerformQuadratureOfObjectiveFunctionSecondTerm(points_set, m_index, a, b):
    
    def DObjectiveFunctionDt(t):
        ais = [points_set.AsTs[m_index][0][i] for i in range(len(points_set.AsTs[m_index][0]))]
        tis = [points_set.AsTs[m_index][1][i] for i in range(len(points_set.AsTs[m_index][0]))]
        alphas = [math.sqrt(math.exp(1) / ti) for ti in tis]
        betas = [- 0.5 / ti for ti in tis]
        return abs(- 0.5 * t ** (- 1.5) - sum([ais[i] * alphas[i] * betas[i] * math.exp(betas[i] * t) for i in range(len(points_set.AsTs[m_index][0]))]))
    
    from scipy.integrate import quad
    value, error = quad(DObjectiveFunctionDt, a, b, limit = 5000)
    return value, error
def ObjectiveFunction(points_set, m_index, pp):
    lower_limit_for_tail = 200000
    return 1.0 / math.sqrt(pp.t_w) * (abs(1 - K(points_set, m_index, 1, pp)) + PerformQuadratureOfObjectiveFunctionSecondTerm(points_set, m_index, 1, lower_limit_for_tail)[0] + ExactErrorTail(points_set, m_index, lower_limit_for_tail))   

class HinsbergPointsSetGivenNorm:
    def __init__(self, code):
        self.code = code
        self.AsTs = []
        self.Errors = []
        self.ErrorBounds = []
        
pp = ProblemParameters()
pp.final_time = 10.0
t_w_min = 0.00001
n_doublings = 6
t_ws = [t_w_min * 10 ** k for k in range(n_doublings)]
t_ws = [1.0]
k_calc = 1    
k_max = len(t_ws)


for t_w in t_ws:
    print('Calculation ' + str(k_calc) + ' of ' + str(k_max) + '.')
    print('t_w = ', t_w, '...')
    print()
    pp.t_w = t_w
    pp.final_time_minus_tw = pp.final_time - pp.t_w
    
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
    size_factor = 3
    line_width = (size_factor * k_calc / k_max) ** 1.2
    small_marker_size = 2 + 2 * (size_factor * k_calc / k_max) ** 1.2
    big_marker_size = 2 + 3 * (size_factor * k_calc / k_max) ** 1.2
    maker_width = 0.75 * (size_factor * k_calc) ** 1.2
    
    if k_calc == min(3, k_max):
        plt.plot(t_norm_set.exponential_numbers, t_norm_set.Errors, markersize = small_marker_size, linewidth = line_width, color='b', marker = 'o', label= 't-norm')
        plt.plot(abs_norm_set.exponential_numbers, abs_norm_set.Errors, markersize = big_marker_size, linewidth = line_width, color='k', marker = '*', label= 'abs-value-norm')
        plt.plot(hinsberg_set.exponential_numbers, hinsberg_set.Errors, markersize = 20, mew = maker_width, linewidth = k_calc, color='g', marker = 'x', label= 'points by van Hinsberg et al.')

        plt.plot(t_norm_set.exponential_numbers, t_norm_set.ErrorBounds, color='b', linestyle='--', label= 't-norm (bound)')
        plt.plot(abs_norm_set.exponential_numbers, abs_norm_set.ErrorBounds, color='k', linestyle=':', label= 'abs-value-norm (bound)')
        plt.plot(hinsberg_set.exponential_numbers, hinsberg_set.ErrorBounds, markersize = 10, color='g', linestyle='--', marker = 'v', label= 'points by van Hinsberg et al. (bound)')
    else:
        plt.plot(t_norm_set.exponential_numbers, t_norm_set.Errors, markersize = small_marker_size, linewidth = line_width, color='b', marker = 'o')
        plt.plot(abs_norm_set.exponential_numbers, abs_norm_set.Errors, markersize = big_marker_size, linewidth = line_width, color='k', marker = '*')
        plt.plot(hinsberg_set.exponential_numbers, hinsberg_set.Errors, markersize = 20, mew = maker_width, linewidth = k_calc, color='g', marker = 'x')

        plt.plot(t_norm_set.exponential_numbers, t_norm_set.ErrorBounds, color='b', linestyle='--')
        plt.plot(abs_norm_set.exponential_numbers, abs_norm_set.ErrorBounds, color='k', linestyle=':')
        plt.plot(hinsberg_set.exponential_numbers, hinsberg_set.ErrorBounds, markersize = 10, color='g', linestyle='--', marker = 'v')    
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
plt.savefig('exact_error_with_sinus_min_t_w=' + str(t_w_min) + '_max_tw=' + str(pp.t_w) + '.eps', format='eps', dpi=1000)
plt.show()
