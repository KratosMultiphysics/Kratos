'''
This file plots data from a file continuously.
E.g. the evolution of the coupling iterations, a points displacement, ... can be observed
Philipp Bucher, 15.10.2017
Chair of Structural Analysis, Technical University of Munich
'''

import matplotlib.pyplot as plt
from numpy import loadtxt

class ResultInfoContainer:
    def __init__(self,
                 ColumnIndex=-1,
                 Label="default_lable",
                 YAxisUpperLimit=None,
                 YAxisLowerLimit=None,
                 Factor=1):
        if ColumnIndex == -1:
            raise Exception("Please specify a column index for " + Label)
        self.Label = Label
        self.ColumnIndex = ColumnIndex
        self.YAxisUpperLimit = YAxisUpperLimit # default is "None" if no limit is set
        self.YAxisLowerLimit = YAxisLowerLimit # default is "None" if no limit is set
        self.Factor = Factor

# =============================================================================
file_name     = "node_output.dat"
ref_file_name = "Mok_Results_ref.dat"
num_rows_to_skip   = 1 # in case there is a header in the file
num_points_to_plot = 500000 # number of data points to plot => e.g. seconds to display / delta_t
plot_update_time   = 10 # [sec]

index_x_axis = 0
label_x_axis  = "Time [sec]"

results = []
results.append(ResultInfoContainer(ColumnIndex=1, Label="Displacement X [m]"))
results.append(ResultInfoContainer(ColumnIndex=2, Label="Displacement Y [m]"))
results.append(ResultInfoContainer(ColumnIndex=3, Label="Displacement Z [m]"))

plot_title   = "WINSENT Beam Deformations Tip"

line_style = 'b-'
# =============================================================================

print("===========================================================")
print("INFO: This will run forever, has to be terminated manually!")
print("===========================================================\n")

# Assembling the results
num_results = len(results)
result_indices = ()
for res in results:
    result_indices += (res.ColumnIndex,)

result_labels = ()
for res in results:
    result_labels += (res.Label,)

result_limits = ()
for res in results:
    result_limits += ([res.YAxisLowerLimit, res.YAxisUpperLimit],)

result_factors = ()
for res in results:
    result_factors += (res.Factor,)

col_tuple = (index_x_axis, ) + result_indices
label_tuple = (label_x_axis, ) + result_labels

plt.ion() # this is responsible for the continuous plot updates
num_points_to_plot = int((-1)*num_points_to_plot) # conversion to make the list slicing work

fig,ax = plt.subplots(num_results,1)

if ref_file_name == "":
    using_ref_file = False
else:
    using_ref_file = True

while(True): # You have to kill this manually!
    try:
        data_results = loadtxt(file_name, skiprows=num_rows_to_skip, usecols=col_tuple, unpack=True)
        if using_ref_file:
            data_res_ref = loadtxt(ref_file_name, skiprows=num_rows_to_skip, usecols=col_tuple, unpack=True)
    except IndexError:
        raise Exception("Loading the results failed, check the requested ColumnIndices!")

    plt.gca().cla() # clear axis to update them
    for res_index in range(1,num_results+1):
        if num_results == 1:
            cur_plot = ax
        else:
            cur_plot = ax[res_index-1]
        cur_plot.clear()
        if result_limits[res_index-1][0] is not None: # set lower axis limit
            cur_plot.set_ylim(bottom=result_limits[res_index-1][0])
        if result_limits[res_index-1][1] is not None: # set upper axis limit
            cur_plot.set_ylim(top=result_limits[res_index-1][1])

        current_num_results = len(data_results[res_index])
        plot_start_index = max(0, current_num_results+num_points_to_plot)
        plot_end_index = current_num_results

        data_results[res_index][plot_start_index:plot_end_index] *= result_factors[res_index-1] # premultiply with factor
        cur_plot.plot(data_results[0][plot_start_index:plot_end_index],data_results[res_index][plot_start_index:plot_end_index], line_style, label='New Result')
        if using_ref_file:
            data_res_ref[res_index][plot_start_index:plot_end_index] *= result_factors[res_index-1] # premultiply with factor
            cur_plot.plot(data_res_ref[0][plot_start_index:plot_end_index],data_res_ref[res_index][plot_start_index:plot_end_index], 'r-', label='Reference Result')

        # Adding labels and Title
        cur_plot.set_ylabel(label_tuple[res_index])
        if res_index == 1:
            cur_plot.set_title(plot_title)
        if res_index == num_results: # add x-lable under the last subplot
            cur_plot.set_xlabel(label_tuple[0])

    if using_ref_file:
        plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
    plt.draw()
    plt.pause(plot_update_time)
