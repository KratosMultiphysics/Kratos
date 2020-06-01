# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

@ExaquteTask(returns=2)
def updateGlobalMonteCarloIndexOrder2Dimension0_Task(globalIndexPowerSum_1,globalIndexPowerSum_2,batchIndexPowerSum_1,batchIndexPowerSum_2):
    globalIndexPowerSum_1 = globalIndexPowerSum_1 + batchIndexPowerSum_1
    globalIndexPowerSum_2 = globalIndexPowerSum_2 + batchIndexPowerSum_2
    return globalIndexPowerSum_1,globalIndexPowerSum_2

@ExaquteTask(returns=5)
def updateGlobalMonteCarloIndexOrder2Dimension1_Task(globalIndexPowerSum_10,globalIndexPowerSum_01,globalIndexPowerSum_20,globalIndexPowerSum_11,globalIndexPowerSum_02,batchIndexPowerSum_10,batchIndexPowerSum_01,batchIndexPowerSum_20,batchIndexPowerSum_11,batchIndexPowerSum_02):
    globalIndexPowerSum_10 = globalIndexPowerSum_10 + batchIndexPowerSum_10
    globalIndexPowerSum_01 = globalIndexPowerSum_01 + batchIndexPowerSum_01
    globalIndexPowerSum_20 = globalIndexPowerSum_20 + batchIndexPowerSum_20
    globalIndexPowerSum_11 = globalIndexPowerSum_11 + batchIndexPowerSum_02
    globalIndexPowerSum_02 = globalIndexPowerSum_02 + batchIndexPowerSum_11
    return globalIndexPowerSum_10,globalIndexPowerSum_01,globalIndexPowerSum_20,globalIndexPowerSum_11,globalIndexPowerSum_02

@ExaquteTask(returns=2)
def updateGlobalMonteCarloIndexTimePowerSumsOrder2Dimension0_Task(globalIndexPowerSum_1,globalIndexPowerSum_2,batchIndexPowerSum_1,batchIndexPowerSum_2):
    globalIndexPowerSum_1 = globalIndexPowerSum_1 + batchIndexPowerSum_1
    globalIndexPowerSum_2 = globalIndexPowerSum_2 + batchIndexPowerSum_2
    return globalIndexPowerSum_1,globalIndexPowerSum_2

@ExaquteTask(returns=4)
def updateGlobalMonteCarloIndexTimePowerSumsOrder2Dimension1_Task(globalIndexPowerSumUpper_1,globalIndexPowerSumLower_1,globalIndexPowerSumUpper_2,globalIndexPowerSumLower_2,batchIndexPowerSumUpper_1,batchIndexPowerSumLower_1,batchIndexPowerSumUpper_2,batchIndexPowerSumLower_2):
    globalIndexPowerSumUpper_1 = globalIndexPowerSumUpper_1 + batchIndexPowerSumUpper_1
    globalIndexPowerSumUpper_2 = globalIndexPowerSumUpper_2 + batchIndexPowerSumUpper_2
    globalIndexPowerSumLower_1 = globalIndexPowerSumLower_1 + batchIndexPowerSumLower_1
    globalIndexPowerSumLower_2 = globalIndexPowerSumLower_2 + batchIndexPowerSumLower_2
    return globalIndexPowerSumUpper_1,globalIndexPowerSumLower_1,globalIndexPowerSumUpper_2,globalIndexPowerSumLower_2