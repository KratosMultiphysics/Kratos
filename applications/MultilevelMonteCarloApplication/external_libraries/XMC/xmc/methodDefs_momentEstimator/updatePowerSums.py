# Import PyCOMPSs
from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
# from exaqute.ExaquteTaskLocal import *      # to execute with python3

from xmc.tools import packedList


def updatePowerSumsOrder2Dimension0(power_sum_1,power_sum_2,*samples):
    samples = packedList(samples)
    if (type(samples) is tuple):
        samples = [samples]
    for i in range (len(samples)):
        sample = samples[i][0]
        if (power_sum_1==None):
            power_sum_1 = sample**(1)
            power_sum_2 = sample**(2)
        else:
            power_sum_1 = power_sum_1 + sample**(1)
            power_sum_2 = power_sum_2 + sample**(2)
    return power_sum_1,power_sum_2

updatePowerSumsOrder2Dimension0_Task = ExaquteTask(returns=2)(updatePowerSumsOrder2Dimension0)

def updatePowerSumsOrder4Dimension0(power_sum_1,power_sum_2,power_sum_3,power_sum_4,*samples):
    samples = packedList(samples)
    if (type(samples) is tuple):
        samples = [samples]
    for i in range (len(samples)):
        sample = samples[i][0]
        if (power_sum_1==None):
            power_sum_1 = sample**(1)
            power_sum_2 = sample**(2)
            power_sum_3 = sample**(3)
            power_sum_4 = sample**(4)
        else:
            power_sum_1 = power_sum_1 + sample**(1)
            power_sum_2 = power_sum_2 + sample**(2)
            power_sum_3 = power_sum_3 + sample**(3)
            power_sum_4 = power_sum_4 + sample**(4)
    return power_sum_1,power_sum_2,power_sum_3,power_sum_4

updatePowerSumsOrder4Dimension0_Task = ExaquteTask(returns=4)(updatePowerSumsOrder4Dimension0)

def updatePowerSumsOrder10Dimension0(power_sum_1,power_sum_2,power_sum_3,power_sum_4,power_sum_5,power_sum_6,power_sum_7,power_sum_8,power_sum_9,power_sum_10,*samples):
    samples = packedList(samples)
    if (type(samples) is tuple):
        samples = [samples]
    for i in range (len(samples)):
        sample = samples[i][0]
        if (power_sum_1==None):
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
    return power_sum_1,power_sum_2,power_sum_3,power_sum_4,power_sum_5,power_sum_6,power_sum_7,power_sum_8,power_sum_9,power_sum_10

updatePowerSumsOrder10Dimension0_Task = ExaquteTask(returns=10)(updatePowerSumsOrder10Dimension0)

def updatePowerSumsOrder2Dimension1(power_sum_10,power_sum_01,power_sum_20,power_sum_11,power_sum_02,*samples):
    samples = packedList(samples)
    if (type(samples) is tuple):
        samples = [samples]
    for i in range (len(samples)):
        sample_sum = samples[i][0]+samples[i][1]
        sample_dif = samples[i][0]-samples[i][1]
        if (power_sum_10==None):
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
    return power_sum_10,power_sum_01,power_sum_20,power_sum_11,power_sum_02

updatePowerSumsOrder2Dimension1_Task = ExaquteTask(returns=5)(updatePowerSumsOrder2Dimension1)