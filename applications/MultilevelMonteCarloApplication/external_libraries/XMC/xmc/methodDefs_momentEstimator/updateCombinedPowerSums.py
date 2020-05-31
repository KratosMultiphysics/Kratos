# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

from xmc.tools import packedList

@ExaquteTask(samples={Type: COLLECTION_IN, Depth: 4},returns=3)
def updatePowerSumsOrder2Dimension0_Task(old_sample_counter,samples,power_sum_1,power_sum_2):
    sample_counter = old_sample_counter
    if (type(samples) is tuple):
        samples = [samples]
    for i in range (len(samples)):
        sample = samples[i][0]
        if (power_sum_1 == None):
            power_sum_1 = sample[0][0]
            power_sum_2 = sample[1][0]
            sample_counter = sample_counter + sample[2]
        else:
            power_sum_1 = power_sum_1 + sample[0][0]
            power_sum_2 = power_sum_2 + sample[1][0]
            sample_counter = sample_counter + sample[2]
    return sample_counter,power_sum_1,power_sum_2

@ExaquteTask(samples={Type: COLLECTION_IN, Depth: 4},returns=5)
def updatePowerSumsOrder2Dimension1_Task(old_sample_counter,samples,power_sum_upper_1,power_sum_lower_1,power_sum_upper_2,power_sum_lower_2):
    sample_counter = old_sample_counter
    if (type(samples) is tuple):
        samples = [samples]
    for i in range (len(samples)):
        sample_upper = samples[i][0]
        if (type(samples[i][1]) is list): # index > 0
            sample_lower = samples[i][1]
        else: # index == 0
            sample_lower = [[0.0],[0.0]]
        if (power_sum_upper_1 == None):
            power_sum_upper_1 = sample_upper[0][0]
            power_sum_upper_2 = sample_upper[1][0]
            power_sum_lower_1 = sample_lower[0][0]
            power_sum_lower_2 = sample_lower[1][0]
            sample_counter = sample_counter + sample_upper[2]
        else:
            power_sum_upper_1 = power_sum_upper_1 + sample_upper[0][0]
            power_sum_upper_2 = power_sum_upper_2 + sample_upper[1][0]
            power_sum_lower_1 = power_sum_lower_1 + sample_lower[0][0]
            power_sum_lower_2 = power_sum_lower_2 + sample_lower[1][0]
            sample_counter = sample_counter + sample_upper[2]
    return sample_counter,power_sum_upper_1,power_sum_lower_1,power_sum_upper_2,power_sum_lower_2