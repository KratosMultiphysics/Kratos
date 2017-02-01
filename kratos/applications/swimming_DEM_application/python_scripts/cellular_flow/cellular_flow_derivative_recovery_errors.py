import math
import matplotlib.pyplot as plt

regular_mesh = False

show_math_deriv_or_laplacian = 'L' # 'M' or 'L'
n_divs = [10, 20]

if regular_mesh:
    sizes = [1.0 / n_div for n_div in n_divs]
else:
    sizes = sizes = [0.05, 0.025]

mat_deriv_recovery_types = [1, 2, 3, 4, 5, 6, 7]
laplacian_recovery_types = [1, 2, 3, 4, 6, 7]

def CalculateLastSlopes(sizes, results):
    Delta_result = math.log(results[-1]/results[-2])
    Delta_size   = math.log(sizes[-1]/sizes[-2])
    slope = abs(Delta_result/Delta_size)
    return slope

def FillVectors(recovery_type, average_errors, max_errors, laplacian_or_mat_deriv):
    if laplacian_or_mat_deriv == 'M':
        recovery_tag = 'mat_deriv_errors_'
    else:
        recovery_tag = 'laplacian_errors_'

    for i in range(len(average_errors)):
        average_errors[i] = 0.0
        max_errors[i] = float('-inf')

    for i_size, (n_div, size) in enumerate(zip(n_divs, sizes)):
        if regular_mesh:
            mesh_tag = 'ndiv_' + str(n_div)
        else:
            mesh_tag = 'h_' + str(size)

        with open('errors_recorded/' + recovery_tag + mesh_tag + '_type_' + str(recovery_type) + '.txt', 'r') as errors_file:
            for i_line, line in enumerate(errors_file):
                numbers_str = line.split()
                error = float(numbers_str[0])
                average_errors[i_size] += error
                max_errors[i_size] = max(max_errors[i_size], error)
                average_errors[i_size] /= i_line + 1

    return average_errors, max_errors



min_error = float("inf")
max_error = - float("inf")

if show_math_deriv_or_laplacian == 'M':
    mat_deriv_average_errors = [0. for size in sizes]
    mat_deriv_max_errors = [float('-inf') for size in sizes]
    for method in mat_deriv_recovery_types:
        FillVectors(method, mat_deriv_average_errors, mat_deriv_max_errors, 'M')
        if method == 1:
            mat_deriv_type = 'standard'
            line_width = 1
            color = 'r'
            marker_type = '-*'
            marker_size = 10

        elif method == 2:
            mat_deriv_type = 'Zhang and Naga (2005)'
            line_width = 1
            color = 'k'
            marker_type = '-o'
            marker_size = 10

        elif method == 3:
            mat_deriv_type = 'L2 (lumped)'
            line_width = 1
            color = 'b'
            marker_type = '->'
            marker_size = 10

        elif method == 4:
            mat_deriv_type = 'L2'
            line_width = 1
            color = 'g'
            marker_type = '->'
            marker_size = 10

        elif method == 5:
            mat_deriv_type = 'L2 only gradient'
            line_width = 1
            color = 'c'
            marker_type = '-v'
            marker_size = 10

        elif method == 6:
            mat_deriv_type = 'Fortin et al. (2012)'
            line_width = 1
            color = 'brown'
            marker_type = '-v'
            marker_size = 10

        elif method == 7:
            mat_deriv_type = 'Zhang and Naga (2005)'
            line_width = 1
            color = 'm'
            marker_type = '-v'
            marker_size = 10
        try:
            mat_deriv_final_slope = CalculateLastSlopes(sizes, mat_deriv_average_errors)
            mat_deriv_slope_msg = ' (m = ' + str(round(mat_deriv_final_slope, 2)) + ')'
        except:
            mat_deriv_slope_msg = ''

        expected_order = 2
        min_error = min(min_error, mat_deriv_average_errors[-1])
        max_error = max(max_error, mat_deriv_average_errors[0])
        plt.plot(sizes, mat_deriv_average_errors, '-v', ms = marker_size, color=color, label= mat_deriv_type + mat_deriv_slope_msg, linewidth = line_width, linestyle='solid', markersize = 20)
        #plt.plot(sizes, mat_deriv_max_errors,'-*', color=color, label= mat_deriv_type + ' material derivative (maximum)', linewidth = 2 * line_width, linestyle='dashed', markersize = 20)

elif show_math_deriv_or_laplacian == 'L':
    laplacian_average_errors = [0. for size in sizes]
    laplacian_max_errors = [float('-inf') for size in sizes]
    for method in laplacian_recovery_types:
        FillVectors(method, laplacian_average_errors, laplacian_max_errors, 'L')
        if method == 1:
            laplacian_type = 'standard'
            line_width = 1
            color = 'r'
            marker_type = '-*'
            marker_size = 10

        elif method == 2:
            laplacian_type = 'Zhang and Naga (2005)'
            line_width = 1
            color = 'k'
            marker_type = '-o'
            marker_size = 10

        elif method == 3:
            laplacian_type = 'L2 divergence of gradient; gradient L2 (lumped)'
            line_width = 1
            color = 'b'
            marker_type = '->'
            marker_size = 10

        elif method == 4:
            laplacian_type = 'L2 divergence of gradient; gradient L2'
            line_width = 1
            color = 'g'
            marker_type = '->'
            marker_size = 10

        elif method == 5:
            laplacian_type = 'L2 divergence of gradient; gradient L2'
            line_width = 1
            color = 'c'
            marker_type = '-v'
            marker_size = 10

        elif method == 6:
            laplacian_type = 'L2 divergence of gradient'
            line_width = 1
            color = 'brown'
            marker_type = '-v'
            marker_size = 10

        elif method == 7:
            laplacian_type = 'Guo (2016)'
            line_width = 1
            color = 'm'
            marker_type = '-v'
            marker_size = 10
        try:
            laplacian_final_slope = CalculateLastSlopes(sizes, laplacian_average_errors)
            laplacian_slope_msg = ' (m = ' + str(round(laplacian_final_slope, 2)) + ')'
        except:
            laplacian_slope_msg = ''

        expected_order = 1
        min_error = min(min_error, laplacian_average_errors[-1])
        max_error = max(max_error, laplacian_average_errors[0])
        plt.plot(sizes, laplacian_average_errors,'-^', color=color, label= laplacian_type + laplacian_slope_msg, linewidth = line_width, linestyle='solid', markersize = 20)
    #plt.plot(sizes, laplacian_max_errors,'-^', color=color, label= laplacian_type + ' laplacian (maximum)', linewidth = 2 * line_width, linestyle='dashed', markersize = 20)
plt.semilogy()
plt.semilogx()
plt.axis('equal')


min_error /= 2

if regular_mesh:
    slope = [min_error * (n_divs[-1] / n_div) ** expected_order for n_div in n_divs]
    plot_name = 'derivative_recovery_errors_regular.pdf'
else:
    slope = [min_error * (size / sizes[-1]) ** expected_order for size in sizes]
    plot_name = 'derivative_recovery_errors_irregular.pdf'
plt.plot(sizes, slope, linestyle='dashed',  label='slope = ' + str(expected_order))
plt.ylim((min_error / 10, max_error * 10))
plt.xlabel('$h$', fontsize=20)

if show_math_deriv_or_laplacian == 'M':
    plt.ylabel('$E_1$', fontsize=20)
else:
    plt.ylabel('$E_2$', fontsize=20)

plt.legend(loc = 'upper left')
plt.savefig(plot_name)
plt.show()
