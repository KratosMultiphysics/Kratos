import math
import matplotlib.pyplot as plt

regular_mesh = True
n_divs = [10, 20, 40, 80]
if regular_mesh:
    sizes = [1.0 / n_div for n_div in n_divs]
else:
    sizes = [0.062, 0.031, 0.016, 0.008]

recovery_types = [1, 2, 3, 4, 5]

def FillVectors(recovery_type, mat_deriv_average_errors, mat_deriv_max_errors, laplacian_average_errors, laplacian_max_errors):
    for i in range(len(mat_deriv_average_errors)):
        mat_deriv_average_errors[i] = 0.0
        mat_deriv_max_errors[i] = float('-inf')
        laplacian_average_errors[i] = 0.0
        laplacian_max_errors[i] = float('-inf')
    i_size = 0
    for size in sizes:
        'laplacian_errors_h_0.025_mat_deriv_type_2_lapl_type_1'
        if regular_mesh:
            mesh_tag = 'ndiv_' + str(n_divs[i_size])
        else:
            mesh_tag = 'h_' + + str(size)
            
        with open('errors_recorded/mat_deriv_errors_' + mesh_tag + '_type_' + str(recovery_type) + '.txt', 'r') as errors_file:
            i_line = 0
            for line in errors_file:                
                i_line += 1
                numbers_str = line.split()
                error = float(numbers_str[0])
                mat_deriv_average_errors[i_size] += error
                mat_deriv_max_errors[i_size] = max(mat_deriv_max_errors[i_size], error)
        mat_deriv_average_errors[i_size] /= i_line

        with open('errors_recorded/laplacian_errors_' + mesh_tag + '_type_' + str(recovery_type) + '.txt', 'r') as errors_file:
            i_line = 0
            for line in errors_file:                
                i_line += 1
                numbers_str = line.split()
                error = float(numbers_str[0])
                laplacian_average_errors[i_size] += error
                laplacian_max_errors[i_size] = max(laplacian_max_errors[i_size], error)
        laplacian_average_errors[i_size] /= i_line
        i_size += 1 
    print('mat_deriv', mat_deriv_average_errors)
        
    return mat_deriv_average_errors, mat_deriv_max_errors, laplacian_average_errors, laplacian_max_errors

mat_deriv_average_errors = [0. for size in sizes]
mat_deriv_max_errors = [float('-inf') for size in sizes]
laplacian_average_errors = [0. for size in sizes]
laplacian_max_errors = [float('-inf') for size in sizes]     
min_error = float("inf")

for method in recovery_types:
    FillVectors(method, mat_deriv_average_errors, mat_deriv_max_errors, laplacian_average_errors, laplacian_max_errors)
    if method == 1:
        mat_deriv_type = 'standard'
        line_width = 1
        color = 'r'
        marker_type = '--'
        marker_size = 20
    elif method == 2:
        mat_deriv_type = 'superconvergent'
        line_width = 1
        color = 'k'
        marker_type = '-+'
        marker_size = 20
        
    elif method == 3:
        mat_deriv_type = 'L2 (lumped)'
        line_width = 1
        color = 'b'
        marker_type = '-*'
        marker_size = 5
        
    elif method == 4:
        mat_deriv_type = 'L2'
        line_width = 1
        color = 'b'
        marker_type = '-*'
        marker_size = 15
        
    elif method == 5:
        mat_deriv_type = 'L2 only gradient'
        line_width = 1
        color = 'b'
        marker_type = '-*'
        marker_size = 30

    plt.plot(sizes, mat_deriv_average_errors, marker_type, ms = marker_size, color=color, label= mat_deriv_type + ' material derivative (average)', linewidth = line_width, linestyle='solid', markersize = 20)
    #plt.plot(sizes, mat_deriv_max_errors, marker_type, ms = marker_size, color=color, label= mat_deriv_type + ' material derivative (maximum)', linewidth = 2 * line_width, linestyle='dashed', markersize = 20)
    #plt.plot(sizes, laplacian_average_errors,'-^', color=color, label= mat_deriv_type + ' laplacian (average)', linewidth = line_width, linestyle='solid', markersize = 20)
    #plt.plot(sizes, laplacian_max_errors,'-^', color=color, label= mat_deriv_type + ' laplacian (maximum)', linewidth = 2 * line_width, linestyle='dashed', markersize = 20)
    min_error = min(min_error, mat_deriv_average_errors[-1], laplacian_average_errors[-1])
    plt.semilogy()
    plt.semilogx()
    plt.axis('equal')
min_error /= 2
slope_1 = [min_error * 10 * (n_divs[-1] / n_div) for n_div in n_divs]
slope_2 = [min_error * (n_divs[-1] / n_div) ** 2 for n_div in n_divs]
plt.plot(sizes, slope_1, linestyle='dashed',  label='slope = 1')
plt.plot(sizes, slope_2, linestyle='dashed',  label='slope = 2')
plt.legend(loc = 'upper left')
plt.savefig('derivative_recovery_errors.pdf')    
plt.show()