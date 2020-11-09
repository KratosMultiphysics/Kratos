import sys
sys.dont_write_bytecode = True

from xmc.distributedEnvironmentFramework import *

import numpy as np

import xmc

@ExaquteTask(returns=1)
def abs_Task(value):
    return abs(value)

@ExaquteTask(returns=1)
def multiplyByScalar_Task(value, scalar):
    return value*scalar

if __name__ == "__main__":

    # RandomGeneratorWrapper
    randomGeneratorInputDict = {'parameters': [0, 1],'generator': 'xmc.methodDefs_randomGeneratorWrapper.generator.normal'}

    # SolverWrapper
    x = [0.11,1.8,0.4,2.1,0.0003,-1.87]
    solverWrapperInputDict = {'index':[],'parameters':x,'outputDimension':1}
    # Resulting distribution parameters:
    # print( x[0]*np.abs((1-np.exp(x[1]))), x[1] )
    # print( (x[2]*(1-np.exp(x[3])))**2, 2*x[3] )
    # print( x[4]*(1-np.exp(x[5])),x[5] )

    # SampleGenerator
    samplerInputDict = {'randomGenerator': 'xmc.randomGeneratorWrapper.RandomGeneratorWrapper',
                'randomGeneratorInputDictionary':randomGeneratorInputDict,
                'solverWrapper': 'xmc.classDefs_solverWrapper.multiLevelRNGSolverWrapper.MultiLevelRNGSolverWrapper',
                'solverWrapperInputDictionary':solverWrapperInputDict,
                'solverWrapperIndices':None}

    # Moment Estimators
    qoiEstimatorInputDict = {'order':1,'indexSetDimension':1,
        "updatedPowerSums":"xmc.methodDefs_momentEstimator.updatePowerSums.updatePowerSumsOrder2Dimension1_Task",
        "centralMomentComputer":"xmc.methodDefs_momentEstimator.computeCentralMoments.centralMomentWrapper_Task",
        "centralMomentErrorComputer":"xmc.methodDefs_momentEstimator.computeErrorEstimation.centralMomentErrorWrapper_Task"}
    costEstimatorInputDict = {'order':1,'indexSetDimension':0,
        "updatedPowerSums":"xmc.methodDefs_momentEstimator.updatePowerSums.updatePowerSumsOrder2Dimension0_Task",
        "centralMomentComputer":"xmc.methodDefs_momentEstimator.computeCentralMoments.centralMomentWrapper_Task",
        "centralMomentErrorComputer":"xmc.methodDefs_momentEstimator.computeErrorEstimation.centralMomentErrorWrapper_Task"}


    # MonteCarloIndex Constructor
    mci_inputDict = {'indexValue': None,
                'sampler':'xmc.sampleGenerator.SampleGenerator',
                'samplerInputDictionary':samplerInputDict,
                'qoiEstimator':['xmc.momentEstimator.MomentEstimator'],
                'qoiEstimatorInputDictionary':[qoiEstimatorInputDict],
                'costEstimator':'xmc.momentEstimator.MomentEstimator',
                'costEstimatorInputDictionary':costEstimatorInputDict,
                'areSamplesRecycled':True
    }

    ################################# RUN TIME GENERATED ENTITIES END HERE #######################

    # MonoCriterion
    refinement_ratio = 2
    final_tol_ratio = 1.3
    initial_tol = 0.5
    final_tol = 0.1

    number_iterations = int(np.ceil(np.log(initial_tol/final_tol)/np.log(refinement_ratio)))
    tolerance_array = [final_tol*(refinement_ratio**(number_iterations-i)) for i in range(number_iterations)]
    tolerance_array.append(final_tol)
    tolerance_array.append(final_tol/final_tol_ratio)
    statErrorCriterion = xmc.monoCriterion.MonoCriterion('xmc.methodDefs_monoCriterion.criterionFunctions.isLowerThanOrEqualTo',tolerance_array)
    exceedLowerBoundCriterion = xmc.monoCriterion.MonoCriterion('xmc.methodDefs_monoCriterion.criterionFunctions.isGreaterThanOrEqualTo',[len(tolerance_array)])
    # TODO Should not criterion function below be 'isLowerThanOrEqualTo'?
    exceedUpperBoundCriterion = xmc.monoCriterion.MonoCriterion('xmc.methodDefs_monoCriterion.criterionFunctions.isGreaterThanOrEqualTo',[20])

    # MultiCriterion
    criteriaDict = [statErrorCriterion, exceedLowerBoundCriterion, exceedUpperBoundCriterion]
    criterionInputs = [['error0'],['algorithmCost'],['algorithmCost']]
    criterion = xmc.multiCriterion.MultiCriterion(criteria=criteriaDict,
                                                  inputsForCriterion=criterionInputs,
                                                  interpreter='xmc.methodDefs_multiCriterion.interpreter.interpretAsConvergenceAndIterationBounds',
                                                  flag='xmc.methodDefs_multiCriterion.flag.plainFlag')


    # ErrorEstimator
    # TODO we should define xmc.methodDefs_errorEstimator.errorEstimation.Variance+BiasError_Task
    MSEErrorEstimator = xmc.errorEstimator.ErrorEstimator(
        error='xmc.methodDefs_errorEstimator.errorEstimation.errorEstimationTerr_Task',
        parameters=[0.95])

    # ModelEstimator
    inputDictModelEstimator = {'valueForParameters':
    'xmc.methodDefs_modelEstimator.valueForParameters.geometricModel',
                  'updater': 'xmc.methodDefs_modelEstimator.update.updatePredictorGeometric_Task'}
    biasPredictor = xmc.modelEstimator.ModelEstimator(**inputDictModelEstimator)
    variancePredictor = xmc.modelEstimator.ModelEstimator(**inputDictModelEstimator)
    costPredictor = xmc.modelEstimator.ModelEstimator(**inputDictModelEstimator)

    # BayesianEstimator
    inputDictBayesianEstimator = {'parameters':[.1,.1],
                    'blendFunction':'xmc.methodDefs_bayesianEstimator.blend.bayesianUpdate'}
    bayesianEstimator = xmc.bayesianEstimator.BayesianEstimator(**inputDictBayesianEstimator)

    # HierarchyOptimiser
    inputDict = {'indexSpace': [10],
                 'toleranceSplittingBounds': [.7]*2,
                 'optimalIndexSet': 'xmc.methodDefs_hierarchyOptimiser.optimalIndexSet.maximalLevelsBiasModel',
                 'optimalSampleNumbers': 'xmc.methodDefs_hierarchyOptimiser.optimalSampleNumbers.multiLevelCostMinimisationTErr',
                 'defaultHierarchy': [[[0],25],[[1],25],[[2],25]],
                 'varianceBlender':bayesianEstimator,
                 'isVarianceBlended':True
                }
    hierarchyCostOptimiser = xmc.hierarchyOptimiser.HierarchyOptimiser(**inputDict)

    # EstimationAssembler
    expectationAssembler = xmc.estimationAssembler.EstimationAssembler(
        assembleEstimation='xmc.methodDefs_estimationAssembler.assembleEstimation.assembleValue_Task')
    biasAssembler = xmc.estimationAssembler.EstimationAssembler(
        assembleEstimation='xmc.methodDefs_estimationAssembler.assembleEstimation.assembleBias_Task')
    varianceAssembler = xmc.estimationAssembler.EstimationAssembler(
        assembleEstimation='xmc.methodDefs_estimationAssembler.assembleEstimation.assembleStatisticalError_Task')

    # MonteCarloSampler
    inputDict = {'indices': [],
                'indexConstructor':'xmc.monteCarloIndex.MonteCarloIndex',
                'indexConstructorDictionary':mci_inputDict,
                'assemblers': [expectationAssembler,biasAssembler,varianceAssembler],
                'estimatorsForAssembler': [ [[0,[1, True, False]]],[[0,[1, True, False]]],[[0,[1, True, True]]] ],
                'qoiPredictor': [biasPredictor,variancePredictor],
                'estimatorsForPredictor': [ [0,[1, True, False],abs_Task] , [0,[1, True, True],multiplyByScalar_Task] ],
                'costEstimatorsForPredictor': [1, True, False],
                'costPredictor': costPredictor,
                'errorEstimators': [MSEErrorEstimator],
                'assemblersForError': [[1,2]],
                'isCostUpdated':True
            }
    mcSampler = xmc.monteCarloSampler.MonteCarloSampler(**inputDict)

    # XMCAlgorithm
    inputDict = {'monteCarloSampler': mcSampler,
                'hierarchyOptimiser': hierarchyCostOptimiser,
                'stoppingCriterion': criterion,
                'errorsForStoppingCriterion': [0],
                'predictorsForHierarchy': [0,1],
                'costPredictorForHierarchy': [1,True,False],
                'estimatorsForHierarchy': [[0,[1,True,False]],[0,[1,True,True]]],
                'costEstimatorForHierarchy': [1,True,False],
                'tolerancesForHierarchy': [0],
                'errorParametersForHierarchy':[0],
        }

    algo = xmc.XMCAlgorithm(**inputDict)

    algo.runXMC()
