from exaqute import *

def updatePowerSumsOrder2Dimension0(samples, power_sum_1, power_sum_2):
    if isinstance(samples,tuple):
        #TODO Probably pointless since deprecation of unpacked lists
        samples = list(samples)
    for i in range(len(samples)):
        sample = samples[i][0]
        if power_sum_1 is None:
            power_sum_1 = sample**(1)
            power_sum_2 = sample**(2)
        else:
            power_sum_1 = power_sum_1 + sample**(1)
            power_sum_2 = power_sum_2 + sample**(2)
    return power_sum_1, power_sum_2

@task(keep=True, returns=2, samples={Type: COLLECTION_IN, Depth: 2})
def updatePowerSumsOrder2Dimension0_Task(samples, *powerSums):
    return updatePowerSumsOrder2Dimension0(samples, *powerSums)

def updatePowerSumsOrder4Dimension0(samples, power_sum_1, power_sum_2, power_sum_3, power_sum_4):
    if isinstance(samples,tuple):
        #TODO Probably pointless since deprecation of unpacked lists
        samples = list(samples)
    for i in range (len(samples)):
        sample = samples[i][0]
        if power_sum_1 is None:
            power_sum_1 = sample**(1)
            power_sum_2 = sample**(2)
            power_sum_3 = sample**(3)
            power_sum_4 = sample**(4)
        else:
            power_sum_1 = power_sum_1 + sample**(1)
            power_sum_2 = power_sum_2 + sample**(2)
            power_sum_3 = power_sum_3 + sample**(3)
            power_sum_4 = power_sum_4 + sample**(4)
    return power_sum_1, power_sum_2, power_sum_3, power_sum_4

@task(keep=True, returns=4, samples={Type: COLLECTION_IN, Depth: 2})
def updatePowerSumsOrder4Dimension0_Task(samples, *powerSums):
    return updatePowerSumsOrder4Dimension0(samples, *powerSums)

def updatePowerSumsOrder10Dimension0(samples, power_sum_1, power_sum_2, power_sum_3, power_sum_4, power_sum_5, power_sum_6, power_sum_7, power_sum_8, power_sum_9, power_sum_10):
    if isinstance(samples,tuple):
        #TODO Probably pointless since deprecation of unpacked lists
        samples = list(samples)
    for i in range (len(samples)):
        sample = samples[i][0]
        if power_sum_1 is None:
            power_sum_1 = sample**(1)
            power_sum_2 = sample**(2)
            power_sum_3 = sample**(3)
            power_sum_4 = sample**(4)
            power_sum_5 = sample**(5)
            power_sum_6 = sample**(6)
            power_sum_7 = sample**(7)
            power_sum_8 = sample**(8)
            power_sum_9 = sample**(9)
            power_sum_10 = sample**(10)
        else:
            power_sum_1 = power_sum_1 + sample**(1)
            power_sum_2 = power_sum_2 + sample**(2)
            power_sum_3 = power_sum_3 + sample**(3)
            power_sum_4 = power_sum_4 + sample**(4)
            power_sum_5 = power_sum_5 + sample**(5)
            power_sum_6 = power_sum_6 + sample**(6)
            power_sum_7 = power_sum_7 + sample**(7)
            power_sum_8 = power_sum_8 + sample**(8)
            power_sum_9 = power_sum_9 + sample**(9)
            power_sum_10 = power_sum_10 + sample**(10)
    return power_sum_1, power_sum_2, power_sum_3, power_sum_4, power_sum_5, power_sum_6, power_sum_7, power_sum_8, power_sum_9, power_sum_10

@task(keep=True, returns=10, samples={Type: COLLECTION_IN, Depth: 2})
def updatePowerSumsOrder10Dimension0_Task(samples, *powerSums):
    return updatePowerSumsOrder10Dimension0(samples, *powerSums)

def updatePowerSumsOrder2Dimension1(samples, power_sum_10, power_sum_01, power_sum_20, power_sum_11, power_sum_02):
    if isinstance(samples,tuple):
        #TODO Probably pointless since deprecation of unpacked lists
        samples = list(samples)
    for i in range (len(samples)):
        sample_sum = samples[i][0]+samples[i][1]
        sample_dif = samples[i][0]-samples[i][1]
        if power_sum_10 is None:
            power_sum_10 = sample_sum**(1)
            power_sum_01 = sample_dif**(1)
            power_sum_20 = sample_sum**(2)
            power_sum_11 = sample_sum**(1) * sample_dif**(1)
            power_sum_02 = sample_dif**(2)
        else:
            power_sum_10 = power_sum_10 + sample_sum**(1)
            power_sum_01 = power_sum_01 + sample_dif**(1)
            power_sum_20 = power_sum_20 + sample_sum**(2)
            power_sum_11 = power_sum_11 + sample_sum**(1) * sample_dif**(1)
            power_sum_02 = power_sum_02 + sample_dif**(2)
    return power_sum_10, power_sum_01, power_sum_20, power_sum_11, power_sum_02

@task(keep=True, returns=5, samples={Type: COLLECTION_IN, Depth: 2})
def updatePowerSumsOrder2Dimension1_Task(samples, *powerSums):
    return updatePowerSumsOrder2Dimension1(samples, *powerSums)
