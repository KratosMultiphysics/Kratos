import math
import matplotlib.pyplot as plt
import h5py

regular_mesh = False
show_math_deriv_or_laplacian = 'M' # 'M' or 'L'
mat_deriv_recovery_types = [1, 3, 6, 7]
laplacian_recovery_types = [1, 3, 6]

marker_size = 10
line_width = 1

def CalculateLastSlopes(sizes, results):
    Delta_result = math.log(results[-1]/results[-2])
    Delta_size   = math.log(sizes[-1]/sizes[-2])
    slope = abs(Delta_result/Delta_size)
    return slope

def FillVectors(recovery_type, sizes, average_errors, max_errors, laplacian_or_mat_deriv):
    del sizes[:]
    del average_errors[:]
    del max_errors[:]
    if laplacian_or_mat_deriv == 'M':
        recovery_tag = 'material derivative/method = ' + str(recovery_type)
    else:
        recovery_tag = 'laplacian/method = ' + str(recovery_type)
    if regular_mesh:
        mesh_tag = 'regular mesh'
    else:
        mesh_tag = 'irregular mesh'
    recovery_tag += '/' + mesh_tag

    with h5py.File('errors_recorded/recovery_errors.hdf5', 'r') as f:
        for i, size_dset in enumerate(f[recovery_tag].values()):
            size = float(size_dset.name.split()[-1])
            if regular_mesh:
                size = 1. / size
            sizes.append(size)
            average_errors.append(size_dset[0])
            max_errors.append(size_dset[1])

min_error = float("inf")
max_error = - float("inf")
average_errors = []
max_errors = []
sizes = []
fig = plt.figure(figsize = (12,10))
ax = fig.add_subplot(1,1,1)

if show_math_deriv_or_laplacian == 'M':
    marker_type = 'v'
    for method in mat_deriv_recovery_types:
        FillVectors(method, sizes, average_errors, max_errors, 'M')
        if method == 1:
            mat_deriv_type = 'standard'
            color = 'r'
        elif method == 2:
            mat_deriv_type = 'Zhang and Naga 2005'
            color = 'k'
        elif method == 3:
            mat_deriv_type = 'L2-lumped'
            color = 'b'
        elif method == 4:
            mat_deriv_type = 'L2'
            color = 'g'
        elif method == 5:
            mat_deriv_type = 'L2 only gradient'
            color = 'c'
        elif method == 6:
            mat_deriv_type = 'Pouliot et al. 2012'
            color = 'brown'
        elif method == 7:
            mat_deriv_type = 'Zhang and Naga 2005'
            color = 'm'
        try:
            mat_deriv_final_slope = CalculateLastSlopes(sizes, average_errors)
            mat_deriv_slope_msg = ' (m = ' + str(round(mat_deriv_final_slope, 2)) + ')'
        except:
            mat_deriv_slope_msg = ''

        expected_order = 2
        min_error = min(min_error, average_errors[-1])
        max_error = max(max_error, average_errors[0])
        plt.plot(sizes, average_errors, marker = marker_type, ms = marker_size, color=color, label= mat_deriv_type + mat_deriv_slope_msg, linewidth = line_width, linestyle='solid', markersize = 20)
        #plt.plot(sizes, mat_deriv_max_errors,'-*', color=color, label= mat_deriv_type + ' material derivative (maximum)', linewidth = 2 * line_width, linestyle='dashed', markersize = 20)

elif show_math_deriv_or_laplacian == 'L':
    marker_type = '^'
    for method in laplacian_recovery_types:
        FillVectors(method, sizes, average_errors, max_errors, 'L')
        if method == 1:
            laplacian_type = 'standard'
            color = 'r'
        elif method == 2:
            laplacian_type = 'Zhang and Naga 2005'
            color = 'k'
        elif method == 3:
            laplacian_type = 'L2 divergence of gradient from L2-lumped'
            color = 'b'
        elif method == 4:
            laplacian_type = 'L2 divergence of gradient from L2'
            color = 'g'
        elif method == 6:
            laplacian_type = 'L2 divergence of gradient from Pouliot et al. 2012'
            color = 'c'
        elif method == 7:
            laplacian_type = 'Guo et al. 2016'
            color = 'm'
        try:
            laplacian_final_slope = CalculateLastSlopes(sizes, average_errors)
            laplacian_slope_msg = ' (m = ' + str(round(laplacian_final_slope, 2)) + ')'
        except:
            laplacian_slope_msg = ''

        expected_order = 1
        min_error = min(min_error, average_errors[-1])
        max_error = max(max_error, average_errors[0])
        ax.plot(sizes, average_errors, marker = marker_type, color=color, label= laplacian_type + laplacian_slope_msg, linewidth = line_width, linestyle='solid', markersize = 20)
    #plt.plot(sizes, laplacian_max_errors,'-^', color=color, label= laplacian_type + ' laplacian (maximum)', linewidth = 2 * line_width, linestyle='dashed', markersize = 20)

plt.semilogx()
plt.semilogy()
plt.axis('equal')
plt.xlim([10 ** -4, 1])

min_error /= 2
if regular_mesh:
    slope = [min_error * (size / sizes[-1]) ** expected_order for size in sizes]
    plot_name = 'derivative_recovery_errors_regular.pdf'
else:
    slope = [min_error * (size / sizes[-1]) ** expected_order for size in sizes]
    plot_name = 'derivative_recovery_errors_irregular.pdf'

ax.plot(sizes, slope, linestyle='dashed',  label='slope = ' + str(expected_order))
plt.ylim((min(slope) / 10, max_error * 10))
plt.xlabel('$h$', fontsize = 20)

if show_math_deriv_or_laplacian == 'M':
    plt.ylabel('$E_1$', fontsize = 20)
else:
    plt.ylabel('$E_2$', fontsize = 20)

plt.legend(loc = 'upper left')
plt.savefig(plot_name, format='eps', bbox_inches = 'tight')
plt.show()
