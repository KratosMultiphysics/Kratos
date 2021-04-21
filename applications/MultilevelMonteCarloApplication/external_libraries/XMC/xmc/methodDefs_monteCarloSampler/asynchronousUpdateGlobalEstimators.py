from collections import Counter
from exaqute import *

@task(keep=True, global_estimators=INOUT,batch_estimators=IN,global_cost_estimator=INOUT,batch_cost_estimator=IN)
def updateGlobalMomentEstimator_Task(global_estimators,batch_estimators,global_cost_estimator,batch_cost_estimator,batch):
    """
    Method updating global estimators with local batch estimators. Power sums and number of realizations are updated.
    """

    number_estimators = len(batch_estimators)

    if batch == 0 or global_cost_estimator.powerSums[0][0] is None:
        # update power sums
        for qoi_index in range (number_estimators):
            if global_estimators[qoi_index].__class__.__name__ in ["MomentEstimator","CombinedMomentEstimator"]:
                for order in range (global_estimators[qoi_index].powerSumsOrder()):
                    for i in range (0,len(global_estimators[qoi_index].powerSums[order])):
                        global_estimators[qoi_index].powerSums[order][i] = batch_estimators[qoi_index].powerSums[order][i]
            elif global_estimators[qoi_index].__class__.__name__ in ["MultiMomentEstimator","MultiCombinedMomentEstimator"]:
                global_estimators[qoi_index]._powerSums = batch_estimators[qoi_index]._powerSums
            else:
                err_msg = "{} is not supported.\n".format(global_estimators[qoi_index].__class__.__name__)
                err_msg += "Available options are: \"MomentEstimator\", \"CombinedMomentEstimator\", \"MultiMomentEstimator\", \"MultiCombinedMomentEstimator\""
                raise Exception(err_msg)
        for order in range (global_cost_estimator.powerSumsOrder()):
            for i in range (0,len(global_cost_estimator.powerSums[order])):
                global_cost_estimator.powerSums[order][i] = batch_cost_estimator.powerSums[order][i]

        # update sample counter
        for qoi_index in range (number_estimators):
            global_estimators[qoi_index]._sampleCounter = batch_estimators[qoi_index]._sampleCounter
        global_cost_estimator._sampleCounter = batch_cost_estimator._sampleCounter

    elif batch > 0:
        # update power sums
        for qoi_index in range (number_estimators):
            if global_estimators[qoi_index].__class__.__name__ in ["MomentEstimator","CombinedMomentEstimator"]:
                for order in range (global_estimators[qoi_index].powerSumsOrder()):
                    for i in range (0,len(global_estimators[qoi_index].powerSums[order])):
                        global_estimators[qoi_index].powerSums[order][i] = global_estimators[qoi_index].powerSums[order][i] + batch_estimators[qoi_index].powerSums[order][i]
            elif (global_estimators[qoi_index].__class__.__name__ in ["MultiMomentEstimator"]) or (global_estimators[qoi_index].__class__.__name__ in ["MultiCombinedMomentEstimator"] and batch_estimators[qoi_index]._indexSetDimension == 0):
                gE = Counter(global_estimators[qoi_index]._powerSums)
                bE = Counter(batch_estimators[qoi_index]._powerSums)
                gE.update(bE)
                global_estimators[qoi_index]._powerSums = gE
            elif global_estimators[qoi_index].__class__.__name__ in ["MultiCombinedMomentEstimator"] and batch_estimators[qoi_index]._indexSetDimension > 0:
                # prepare dictionaries for updating and update them
                gEU = Counter(global_estimators[qoi_index]._powerSums["upper"])
                bEU = Counter(batch_estimators[qoi_index]._powerSums["upper"])
                gEL = Counter(global_estimators[qoi_index]._powerSums["lower"])
                bEL = Counter(batch_estimators[qoi_index]._powerSums["lower"])
                gEU.update(bEU)
                gEL.update(bEL)
                # clear global estimator dictionaries and replace with updated dictionaries
                global_estimators[qoi_index]._powerSums["upper"].clear()
                global_estimators[qoi_index]._powerSums["lower"].clear()
                global_estimators[qoi_index]._powerSums["upper"].update(gEU)
                global_estimators[qoi_index]._powerSums["lower"].update(gEL)
            else:
                err_msg = "{} is not supported.\n".format(global_estimators[qoi_index].__class__.__name__)
                err_msg += "Available options are: \"MomentEstimator\", \"CombinedMomentEstimator\", \"MultiMomentEstimator\", \"MultiCombinedMomentEstimator\""
                raise Exception(err_msg)
        for order in range (global_cost_estimator.powerSumsOrder()):
            for i in range (0,len(global_cost_estimator.powerSums[order])):
                global_cost_estimator.powerSums[order][i] = global_cost_estimator.powerSums[order][i] + batch_cost_estimator.powerSums[order][i]

        # update sample counter
        for qoi_index in range (number_estimators):
            global_estimators[qoi_index]._sampleCounter = global_estimators[qoi_index]._sampleCounter + batch_estimators[qoi_index]._sampleCounter
        global_cost_estimator._sampleCounter = global_cost_estimator._sampleCounter + batch_cost_estimator._sampleCounter
