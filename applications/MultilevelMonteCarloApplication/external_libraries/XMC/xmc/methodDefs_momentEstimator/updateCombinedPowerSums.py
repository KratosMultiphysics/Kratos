from xmc.distributedEnvironmentFramework import *

def updatePowerSumsOrder1Dimension0():
    pass

@ExaquteTask(samples={Type: COLLECTION_IN, Depth: 4},returns=1)
def updatePowerSumsOrder1Dimension0_Task(samples,*args):
    return updatePowerSumsOrder1Dimension0(samples,*args)

def updatePowerSumsOrder2Dimension0(old_sample_counter,samples,power_sum_1,power_sum_2):
    sample_counter = old_sample_counter
    if (type(samples) is tuple):
        samples = [samples]
    for i in range (len(samples)):
        sample = samples[i][0]
        if (power_sum_1 == None):
            power_sum_1 = sample[0][0]
            power_sum_2 = sample[1][0]
            sample_counter = sample_counter + sample[-1]
        else:
            power_sum_1 = power_sum_1 + sample[0][0]
            power_sum_2 = power_sum_2 + sample[1][0]
            sample_counter = sample_counter + sample[-1]
    return sample_counter,power_sum_1,power_sum_2

@ExaquteTask(samples={Type: COLLECTION_IN, Depth: 4},returns=3)
def updatePowerSumsOrder2Dimension0_Task(counter,samples,*args):
    return updatePowerSumsOrder2Dimension0(counter,samples,*args)

def updatePowerSumsOrder10Dimension0(old_sample_counter,samples,power_sum_1,power_sum_2,power_sum_3,power_sum_4,power_sum_5,power_sum_6,power_sum_7,power_sum_8,power_sum_9,power_sum_10):
    sample_counter = old_sample_counter
    if (type(samples) is tuple):
        samples = [samples]
    for i in range (len(samples)):
        sample = samples[i][0]
        if (power_sum_1 == None):
            power_sum_1 = sample[0][0]
            power_sum_2 = sample[1][0]
            power_sum_3 = sample[2][0]
            power_sum_4 = sample[3][0]
            power_sum_5 = sample[4][0]
            power_sum_6 = sample[5][0]
            power_sum_7 = sample[6][0]
            power_sum_8 = sample[7][0]
            power_sum_9 = sample[8][0]
            power_sum_10 = sample[9][0]
            sample_counter = sample_counter + sample[-1]
        else:
            power_sum_1 = power_sum_1 + sample[0][0]
            power_sum_2 = power_sum_2 + sample[1][0]
            power_sum_3 = power_sum_3 + sample[2][0]
            power_sum_4 = power_sum_4 + sample[3][0]
            power_sum_5 = power_sum_5 + sample[4][0]
            power_sum_6 = power_sum_6 + sample[5][0]
            power_sum_7 = power_sum_7 + sample[6][0]
            power_sum_8 = power_sum_8 + sample[7][0]
            power_sum_9 = power_sum_9 + sample[8][0]
            power_sum_10 = power_sum_10 + sample[9][0]
            sample_counter = sample_counter + sample[-1]
    return sample_counter,power_sum_1,power_sum_2,power_sum_3,power_sum_4,power_sum_5,power_sum_6,power_sum_7,power_sum_8,power_sum_9,power_sum_10

@ExaquteTask(samples={Type: COLLECTION_IN, Depth: 4},returns=11)
def updatePowerSumsOrder10Dimension0_Task(counter,samples,*args):
    return updatePowerSumsOrder10Dimension0(counter,samples,*args)

def updatePowerSumsOrder2Dimension1(old_sample_counter,samples,power_sum_upper_1,power_sum_lower_1,power_sum_upper_2,power_sum_lower_2):
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
            sample_counter = sample_counter + sample_upper[-1]
        else:
            power_sum_upper_1 = power_sum_upper_1 + sample_upper[0][0]
            power_sum_upper_2 = power_sum_upper_2 + sample_upper[1][0]
            power_sum_lower_1 = power_sum_lower_1 + sample_lower[0][0]
            power_sum_lower_2 = power_sum_lower_2 + sample_lower[1][0]
            sample_counter = sample_counter + sample_upper[-1]
    return sample_counter,power_sum_upper_1,power_sum_lower_1,power_sum_upper_2,power_sum_lower_2

@ExaquteTask(samples={Type: COLLECTION_IN, Depth: 4},returns=5)
def updatePowerSumsOrder2Dimension1_Task(counter,samples,*args):
    return updatePowerSumsOrder2Dimension1(counter,samples,*args)
