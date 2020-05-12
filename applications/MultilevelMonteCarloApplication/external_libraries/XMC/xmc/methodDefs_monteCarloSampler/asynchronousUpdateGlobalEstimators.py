# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

@ExaquteTask(global_estimators=INOUT,batch_estimators=IN,global_cost_estimator=INOUT,batch_cost_estimator=IN)
def updateGlobalMomentEstimator_Task(global_estimators,batch_estimators,global_cost_estimator,batch_cost_estimator,batch):
    number_estimators = len(batch_estimators)

    if (batch == 0):
        # update power sums
        for qoi_index in range (number_estimators):
            for order in range (global_estimators[qoi_index].powerSumsOrder()):
                for i in range (0,len(global_estimators[qoi_index].powerSums[order])):
                    global_estimators[qoi_index].powerSums[order][i] = batch_estimators[qoi_index].powerSums[order][i]
        for order in range (global_cost_estimator.powerSumsOrder()):
            for i in range (0,len(global_cost_estimator.powerSums[order])):
                global_cost_estimator.powerSums[order][i] = batch_cost_estimator.powerSums[order][i]

        # update sample counter
        for qoi_index in range (number_estimators):
            global_estimators[qoi_index]._sampleCounter = batch_estimators[qoi_index]._sampleCounter
        global_cost_estimator._sampleCounter = batch_cost_estimator._sampleCounter

    elif (batch > 0):
        # update power sums
        for qoi_index in range (number_estimators):
            for order in range (global_estimators[qoi_index].powerSumsOrder()):
                for i in range (0,len(global_estimators[qoi_index].powerSums[order])):
                    global_estimators[qoi_index].powerSums[order][i] = global_estimators[qoi_index].powerSums[order][i] + batch_estimators[qoi_index].powerSums[order][i]
        for order in range (global_cost_estimator.powerSumsOrder()):
            for i in range (0,len(global_cost_estimator.powerSums[order])):
                global_cost_estimator.powerSums[order][i] = global_cost_estimator.powerSums[order][i] + batch_cost_estimator.powerSums[order][i]

        # update sample counter
        for qoi_index in range (number_estimators):
            global_estimators[qoi_index]._sampleCounter = global_estimators[qoi_index]._sampleCounter + batch_estimators[qoi_index]._sampleCounter
        global_cost_estimator._sampleCounter = global_cost_estimator._sampleCounter + batch_cost_estimator._sampleCounter
