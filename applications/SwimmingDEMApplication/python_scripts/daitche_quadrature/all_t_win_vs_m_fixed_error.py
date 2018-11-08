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

def CalculateErrors(points_set, pp, samples = 10):
    Error = 0
    error_bound = ObjectiveFunction(points_set, m, pp)
    tis = [pp.t_w * ti for ti in points_set.AsTs[m][1]]

    for m in points_set.exponential_indices:
        for i in range(n_samples):
            final_time = pp.final_time - 2 * math.pi * i / n_samples
            final_time_minus_tw = final_time - pp.t_w
            ExactIntegral = quad.ExactIntegrationOfSinus(final_time)
            F_w = quad.ExactIntegrationOfSinus(final_time, final_time_minus_tw)
            tis = [pp.t_w * ti for ti in points_set.AsTs[m][1]]
            F_tail = quad.ExactIntegrationOfTail(final_time, final_time_minus_tw, 0, points_set.AsTs[m][0], tis)
            ApproxIntegral = float(F_w) + float(F_tail)
            error_bound = ObjectiveFunction(points_set, m, pp)
            Error += abs((ExactIntegral - ApproxIntegral))
        points_set.Errors.append(Error / samples)
        points_set.ErrorBounds.append(error_bound)

def CalculateError(t_w):
    ExactIntegral = quad.ExactIntegrationOfSinus(pp.final_time)
    F_w = quad.ExactIntegrationOfSinus(pp.final_time, pp.final_time_minus_tw)
    tis = [t_w * ti for ti in points_set.AsTs[m][1]]
    F_tail = quad.ExactIntegrationOfTail(pp.final_time, pp.final_time_minus_tw, 0, points_set.AsTs[m][0], tis)
    ApproxIntegral = float(F_w) + float(F_tail)
    error_bound = ObjectiveFunction(points_set, m, pp)
    Error = abs((ExactIntegral - ApproxIntegral))
    return error, error_bound

def CalculateApproxTailContribution(points_set, m_index, lower_limit):
    if m == 0:
        return 0.
    else:
        ais = points_set.AsTs[m_index][0][:]
        tis = points_set.AsTs[m_index][1][:]
        alphas = [math.sqrt(math.exp(1) / ti) for ti in tis]
        betas = [- 0.5 / ti for ti in tis]
        return sum([ais[i] * alphas[i] * math.exp(betas[i] * lower_limit) for i in range(len(points_set.AsTs[m_index][0]))])

def ExactErrorTail(points_set, m_index, lower_limit):
    ais = points_set.AsTs[m_index][0][:]
    tis = points_set.AsTs[m_index][1][:]
    alphas = [math.sqrt(math.exp(1) / ti) for ti in tis]
    betas = [- 0.5 / ti for ti in tis]
    return 1. / math.sqrt(lower_limit) - sum([ais[i] * alphas[i] * math.exp(betas[i] * lower_limit) for i in range(len(points_set.AsTs[m_index][0]))])

def K(points_set, m_index, t, pp):
    ais = points_set.AsTs[m_index][0][:]
    tis = points_set.AsTs[m_index][1][:]
    alphas = [math.sqrt(math.exp(1) / ti) for ti in tis]
    betas = [- 0.5 / ti for ti in tis]
    return sum([ais[i] * alphas[i] * math.exp(betas[i] * t) for i in range(len(points_set.AsTs[m_index][0]))])

def PerformQuadratureOfObjectiveFunctionSecondTerm(points_set, m_index, a, b):

    def DObjectiveFunctionDt(t):
        ais = points_set.AsTs[m_index][0][:]
        tis = points_set.AsTs[m_index][1][:]
        alphas = [math.sqrt(math.exp(1) / ti) for ti in tis]
        betas = [- 0.5 / ti for ti in tis]
        return abs(- 0.5 * t ** (- 1.5) - sum([ais[i] * alphas[i] * betas[i] * math.exp(betas[i] * t) for i in range(len(points_set.AsTs[m_index][0]))]))

    from scipy.integrate import quad
    value, error = quad(DObjectiveFunctionDt, a, b, limit = 5000)
    return value, error

def ObjectiveFunction(points_set, m_index, pp):
    lower_limit_for_tail = 200000
    return 1.0 / math.sqrt(pp.t_w) * (abs(1 - K(points_set, m_index, 1, pp)) + PerformQuadratureOfObjectiveFunctionSecondTerm(points_set, m_index, 1, lower_limit_for_tail)[0] + ExactErrorTail(points_set, m_index, lower_limit_for_tail))

def sign(x):
    if x > 0:
        return 1
    elif x < 0:
        return -1
    else:
        return 0


def samesign(a, b):
        return a * b > 0

def bisect(func, low, high):
    'Find root of continuous function where f(low) and f(high) have opposite signs'
    f_low = func(low)
    f_high = func(high)
    while samesign(f_low, f_high):
        if f_low < 0:
            high /= 2.
            low /= 2.
        else:
            high *= 2
            low /= 1.2
        f_low = func(low)
        f_high = func(high)

    for i in range(50):
        midpoint = 0.5 * (low + high)
        if samesign(func(low), func(midpoint)):
            low = midpoint
        else:
            high = midpoint
        if (high - low) / high < 1e-3:
            break
        print('tw ', midpoint)

    return midpoint, low, high

def regula_falsi(func, a, b, max_steps=100, tolerance=1e-8):
    p = 0.0
    for loopCount in range(max_steps):
        p = b - (func(b) * ((a-b)/(func(a)-func(b))))
        if math.copysign(func(a), func(b)) != func(a):
            b = p
        else:
            a = p
        if abs(func(p)) < tolerance:
            return p
    print('Root find cancelled at %.9f' % p)

class HinsbergPointsSetGivenNorm:
    def __init__(self, code):
        self.code = code
        self.AsTs = []
        self.Errors = []
        self.ErrorBounds = []

pp = ProblemParameters()
pp.final_time = 2 * math.pi
pp.initial_time = 'MinusInfinity'
ref_errors = [1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4]
# ref_errors = [1e-1, 1e-2, 1e-3]
max_t_w_t_norms = []
max_t_w_t_norm = [10 ** 3] * 11
max_t_w_abs_norm = [10 ** 3] * 11
max_t_w_hinsberg = [10 ** 3]

def CalculateTws(m, max_t_w, points_set, ref_error, n_samples = 40):
    if points_set.code == 'hinsberg_norm':
        m_index = 0
    else:
        m_index = m - 1

    def ObjectiveFunctionMaxTwin(t_w):
        error = 0.
        for i in range(n_samples):
            final_time = 2 * math.pi * (i + 1) / n_samples
            final_time_minus_tw = final_time - t_w
            approximate_tail_contribution = 0.

            if m_index + 1:
                ais = points_set.AsTs[m_index][0][:]
                tis = [t_w * ti for ti in points_set.AsTs[m_index][1]]
                approximate_tail_contribution = float(quad.ExactIntegrationOfTail(final_time = final_time, final_time_minus_tw = final_time_minus_tw, initial_time = 'MinusInfinity', ais = ais, tis = tis))
            exact_tail = quad.ExactIntegrationOfSinus(final_time, a = 'MinusInfinity', b = final_time) - quad.ExactIntegrationOfSinus(final_time, a = final_time_minus_tw, b = final_time)
            error += abs(exact_tail - approximate_tail_contribution)
            # print('a', final_time_minus_tw)
            # print('EXACT', exact_tail)
            # print('APROX',approximate_tail_contribution)
        return error / n_samples / ref_error - 1

    if points_set.code == 'hinsberg_norm':
        m = 0

    high = 10. ** 3
    low = 10. ** - 8
    if points_set.code == 'hinsberg_norm' or m == 0:
        high = min(high, max_t_w[0])
    else:
        high = min(high, 4 * max_t_w[m], 4 * max_t_w[m - 1])

    print(max_t_w)
    max_t_w[m], a, b = bisect(func = ObjectiveFunctionMaxTwin, low = low, high = high)
    # print()
    # print('REGULA!')
    # print()
    # max_t_w[m] = regula_falsi(ObjectiveFunctionMaxTwin, a, b)

t_norm_set = HinsbergPointsSetGivenNorm('t_norm')
abs_norm_set = HinsbergPointsSetGivenNorm('abs_norm')
hinsberg_set = HinsbergPointsSetGivenNorm('hinsberg_norm')

op.FillPoints(t_norm_set, t_norm_set.code)
op.FillPoints(abs_norm_set, abs_norm_set.code)
op.FillPoints(hinsberg_set, hinsberg_set.code)
plt.figure(figsize = (14, 11))

for i_plot, ref_error in enumerate(ref_errors):
    for m in range(0, 11):
        print('*****************************************************************')
        print('m = ', m)
        print('*****************************************************************')
        print()
        print('t-norm')
        print()
        CalculateTws(m, max_t_w_t_norm, t_norm_set, ref_error)
        print()
        print('abs-norm')
        print()
        max_t_w_abs_norm[m] = max_t_w_t_norm[m]
        if m:
            CalculateTws(m, max_t_w_abs_norm, abs_norm_set, ref_error)
        if m == 10:
            print()
            print('hinsberg points')
            print()
            max_t_w_hinsberg[0] = max_t_w_abs_norm[10]
            CalculateTws(m, max_t_w_hinsberg, hinsberg_set, ref_error)

    print('tw tnorm', max_t_w_t_norm)
    print('tw absnorm', max_t_w_abs_norm)
    print('tw hinsberg', max_t_w_hinsberg)
    tw_max = 1.0 #max_t_w_abs_norm[0]
    max_t_w_t_norm = [tw / tw_max for tw in max_t_w_t_norm]
    max_t_w_abs_norm = [tw / tw_max for tw in max_t_w_abs_norm]
    max_t_w_hinsberg = [tw / tw_max for tw in max_t_w_hinsberg]
    # max_t_w_t_norms.append([tw_max, max_t_w_t_norm, max_t_w_abs_norm, max_t_w_hinsberg])
    plt.plot(list(range(11)), max_t_w_t_norm, '-o', color = 'b', ms = 10 * 1.2 ** i_plot, linewidth = 2 * 1.3 ** i_plot, label = '$I_{2t}$')
    plt.plot(list(range(11)), max_t_w_abs_norm, '-*', color = 'k', ms = 13 * 1.2 ** i_plot, linewidth = 2 * 1.3 ** i_plot, label = '$I_{1}$')
    plt.plot([10], max_t_w_hinsberg, 'x', color = 'brown', mew = 5, ms = 10 * 1.2 ** i_plot, linewidth = 2 * 1.3 ** i_plot, label = '$I_{2tH}$')
    plt.semilogy()
    # plt.legend(loc = 'upper right', prop={'size':22.5},frameon=False)
    plt.xlabel('$m$', fontsize = '30')
    percentage = '%g' % (100 * ref_error)
    plt.ylabel('$t_{{w}}$'.format(p = percentage), fontsize = '30')
    ax = plt.gca()
    ax.tick_params(axis='x', pad=20, labelsize=20)
    ax.tick_params(axis='y', pad=10, labelsize=20)
    figure = plt.gcf() # get current figure
    figure.set_size_inches(14, 11)
    plt.xlim([0, 11])
    plt.ylim([0, 10**8])
    plt.tight_layout()
    base_t_w = '%.2f' % tw_max

    # file_name = 'non_dimensional_tw_with_base_tw=' + base_t_w + '_for_' + percentage + '_percent_error.pdf'
    # plt.savefig(file_name, format='pdf', dpi=1000)
    # plt.close()
file_name = 'all_tw_for_x_percent_error.pdf'
#lgnd = plt.legend(loc = 'upper right', prop={'size':22}, frameon=False)
# for handel in lgnd.legendHandles:
#     handel._legmarker.set_markersize(12)
plt.savefig(file_name, format='pdf', dpi=1000)
plt.show()
