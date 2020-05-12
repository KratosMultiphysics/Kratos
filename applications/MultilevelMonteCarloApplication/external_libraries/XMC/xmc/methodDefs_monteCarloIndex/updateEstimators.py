import xmc.tools as tools

# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3
from pycompss.api.parameter import *


@ExaquteTask(qoi_estimators=INOUT,new_samples={Type: COLLECTION_IN, Depth: 4})
def updatePartialQoiEstimators_Task(new_samples,qoi_estimators):
    number_estimators = len(qoi_estimators)
    samples = [[] for _ in range (number_estimators)]
    for sample_aux in new_samples:
        # sample is a list with shape [[QoI_1^l,...,QoI_N^l],[QoI_1^(l-1),...,QoI_N^(l-1)],...]
        # need to convert to shape [[QoI_1^l, QoI_1^(l-1), ...],...,[QoI_N^l, QoI_N^(l-1), ...]]

        # flat the list of lists
        sample=[[val] for lev in sample_aux for val in lev ]
        # if more than one level (MLMC)
        if (len(sample) > number_estimators):
            sample_sublists = tools.splitList(sample,num_sublists=int(len(sample) / number_estimators)) # split list into finer and coarser lists
            sample = [ [] for _ in range (0,number_estimators)]
            for smpl_sublist in sample_sublists:
                for output_counter in range (0,number_estimators):
                    sample[output_counter].append(smpl_sublist[output_counter][0])
        # append to generate [[QoI_1^l, QoI_1^(l-1), ...],...,[QoI_N^l, QoI_N^(l-1), ...]]
        for output_counter in range (0,number_estimators):
            samples[output_counter].append(sample[output_counter])

    # Pass new samples to relevant estimator
    for output_counter in range(number_estimators):
        # Update estimators for a single solver output
        qoi_estimators[output_counter].update(samples[output_counter])

@ExaquteTask(qoi_group=IN, qoi_estimators=INOUT)
def mergeQoiEstimators_Task(qoi_group, qoi_estimators, group_number, group_size):
    number_estimators = len(qoi_estimators)
    for i in range (0, group_size):
        estim_counter = (group_number*group_size)+i
        if estim_counter >= number_estimators:
            break
        else:
            qoi_estimators[estim_counter]=qoi_group[i]

@ExaquteTask(cost_estimator=INOUT,new_times={Type: COLLECTION_IN, Depth: 4})
def updateCostEstimator_Task(new_times,cost_estimator):
    # Pass new times to relevant estimator
    times = []
    for tt in new_times: # loop over samples
        t = tools.summation(*tools.unpackedList(tt))
        times.append([t])
    cost_estimator.update(times)