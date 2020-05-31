import pickle
import sys
import pathlib as pl
import os

# XMC imports
from xmc.tools import dynamicImport
from xmc.tools import instantiateObject
from xmc.tools import splitOneListIntoTwo

# TODO: remove PyCOMPSs import from here
# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

class XMCAlgorithm():
    """
    This top-level class handles the overall algorithm: initialisation as well
    as everything related to error, tolerance and iterations. It also possesses
    the necessary methods and attributes to create new types of algorithms. However,
    the export of results is to be handled outside. Although self-sufficient, this class may inherit from some metaclass.
    """

    def __init__(self, **keywordArgs):
        #TODO - think of a better way to do data dumps and output
        self.isDataDumped = keywordArgs.get('isDataDumped',False)
        self.outputFolderPath = keywordArgs.get('outputFolderPath',None)

        # Attributes
        self.monteCarloSampler = keywordArgs.get('monteCarloSampler')
        self.hierarchyOptimiser = keywordArgs.get('hierarchyOptimiser')
        self.stoppingCriterion = keywordArgs.get('stoppingCriterion')
        self.errorsForStoppingCriterion = keywordArgs.get('errorForStoppingCriterion',None)
        self.predictorsForHierarchy = keywordArgs.get('predictorsForHierarchy',None)
        self.estimatorsForHierarchy = keywordArgs.get('estimatorsForHierarchy',None)
        self.tolerancesForHierarchy = keywordArgs.get('tolerancesForHierarchy',None)
        self.errorParametersForHierarchy = keywordArgs.get('errorParametersForHierarchy',None)
        self.costPredictorForHierarchy = keywordArgs.get('costPredictorForHierarchy',None)
        self.costEstimatorForHierarchy = keywordArgs.get('costEstimatorForHierarchy',None)
        self.positionMaxNumberIterationsCriterion = keywordArgs.get('positionMaxNumberIterationsCriterion',None)
        self.iterationCounter = 0

        # Methods
        self.checkInitialisation = dynamicImport(keywordArgs.get('checkInitialisation','xmc.tools.doNothing'))

    def updateTolerance(self, criteriaToUpdate=None):
        self.stoppingCriterion.updateTolerance(criteriaToUpdate)

    def indexEstimation(self, coordinate, valueMethodArgs):
        return self.monteCarloSampler.indexEstimation(coordinate, valueMethodArgs)

    def indexCostEstimation(self, valueMethodArgs):
        return self.monteCarloSampler.indexCostEstimation(valueMethodArgs)

    def qoiPrediction(self):
        return self.monteCarloSampler.qoiPredictor

    def costPrediction(self):
        return self.monteCarloSampler.costPredictor

    def tolerances(self,*args):
        return self.stoppingCriterion.tolerances(*args)

    def hierarchy(self):
        """
        Returns current hierarchy of the MC estimator.
        """
        return self.monteCarloSampler.hierarchy()

    def splitTolerance(self, splittingParameter):
        """
        Method that interfaces with the MultiCriterion class to apply tolerance
        splitting
        """
        self.stoppingCriterion.splitTolerance(splittingParameter)

    def optimalHierarchy(self):
        """
        Method that interfaces with the HierarchyOptimiser class to compute
        the optimal hierarchy
        """
        old_hierarchy = self.hierarchy()
        input_dict = self.hierarchyOptimiser.inputDictionaryTemplate()
        if(self.iterationCounter == 0):
            newHierarchy = self.hierarchyOptimiser.defaultHierarchy
        else:
            # Build model_list and parameters_list
            model_list = []
            parameter_list = []
            old_parameter_list = []
            qoi_predictors = self.qoiPrediction()
            for qoi_predictor_coordinate in self.predictorsForHierarchy:
                model_list.append(qoi_predictors[qoi_predictor_coordinate]._valueForParameters)
                parameter_list.append(qoi_predictors[qoi_predictor_coordinate].parameters)
                old_parameter_list.append(qoi_predictors[qoi_predictor_coordinate].oldParameters)
            cost_predictor = self.costPrediction()

            # Build estimations_list
            estimation_list = []
            for qoi_estimator_coordinate_and_value_args in self.estimatorsForHierarchy:
                coordinate = qoi_estimator_coordinate_and_value_args[0]
                value_args = qoi_estimator_coordinate_and_value_args[1]
                estimation_list.append(self.indexEstimation(coordinate,value_args))
            if self.costEstimatorForHierarchy is not None:
                cost_estimations = self.indexCostEstimation(self.costEstimatorForHierarchy)
            else:
                cost_estimations = None

            # Build error_parameters_list
            error_parameters_list = []
            # TODO - Triple dereference below!! Add method to get errorEstimator parameters
            # or errorEstimator objects themselves from monteCarloSampler
            if self.errorParametersForHierarchy is not None:
                for coord in self.errorParametersForHierarchy:
                    error_parameters_list.append(self.monteCarloSampler.errorEstimators[coord].parameters)

            input_dict['errorParameters'] = error_parameters_list
            input_dict['oldHierarchy'] = old_hierarchy
            # TODO - Accept these hard coded values from a dictionary somewhere
            input_dict['newSampleNumber'] = 25
            input_dict['tolerances'] = self.tolerances(self.tolerancesForHierarchy)
            input_dict['models'] = model_list
            input_dict['parametersForModel'] = get_value_from_remote(parameter_list)
            input_dict['parametersForModelOld'] = get_value_from_remote(old_parameter_list)
            input_dict['estimations'] = get_value_from_remote(estimation_list)
            if cost_predictor is None:
                input_dict['costModel'] = None
                input_dict['costParameters'] = None
                input_dict['costParametersOld'] = None
            else:
                input_dict['costModel'] = cost_predictor._valueForParameters
                input_dict['costParameters'] = get_value_from_remote(cost_predictor.parameters)
                input_dict['costParametersOld'] = get_value_from_remote(cost_predictor.oldParameters)
            input_dict['costEstimations'] = get_value_from_remote(cost_estimations)
            input_dict['defaultHierarchy'] = self.hierarchyOptimiser.defaultHierarchy

            newHierarchy = self.hierarchyOptimiser.optimalHierarchy(input_dict)
        return newHierarchy,input_dict.get('splittingParameter',None)

    def updateHierarchy(self,newHierarchy):
        """
        Method that interfaces with the monteCarloSample class to execute
        a given hierarchy
        """
        # TODO could be confused with optimalHierarchy. Rename updateSampler or updateSamplerHierarchy?
        self.monteCarloSampler.update(newHierarchy)

    def estimation(self, assemblerCoordinates=None):
        """
        Method that calls the estimation method of monteCarloSampler
        """
        return self.monteCarloSampler.estimation(assemblerCoordinates)

    def errorEstimation(self, errorEstimatorCoordinates=None):
        """
        Method that calls the errorEstimation method of monteCarloSampler
        """
        return self.monteCarloSampler.errorEstimation(errorEstimatorCoordinates)

    def updateHierarchySpace(self,*args):
        """
        Method that interfaces with the HierarchyOptimiser class to compute
        the hierarchy space in which to search for the optimal hierarchy
        """
        self.hierarchyOptimiser.updateHierarchySpace(args)

    def stoppingCriterionFlag(self,currentCost=None):
        """
        Call stoppingCriterion.flag with the proper arguments and return its output
        (a.k.a flag).
        Input argument: currentCost is an indication of the cost the algorithm has entailed
        so far; we usually use the number of iterations.
        Output argument: criterion flag structure as define in the MultiCriterion class.
        """
        # Get errors required for stopping criterion
        errors = get_value_from_remote(self.errorEstimation(self.errorsForStoppingCriterion))

        # Set up dictionary required for stoppingCriterion.flag
        input_dictionary = {}
        for i in range(len(errors)):
            input_dictionary['error'+str(i)] = errors[i]
        input_dictionary['hierarchy'] = self.hierarchy()
        input_dictionary['algorithmCost'] = currentCost

        # Comput flag from dictionary and return
        flag = self.stoppingCriterion.flag(input_dictionary)
        return flag

    def runXMC(self):
        """
        Run an algorithm with generic structure, as described in the documentation.
        Other methods may be implemented in the future, for specific needs that do
        not fit even this generic algorithm.
        """

        self.checkInitialisation()

        # Iteration Loop will start here
        flag = self.stoppingCriterion.flagStructure()
        self.iterationCounter = 0
        #print("Beginning Iterations for tolerance ", self.tolerances(self.errorsForStoppingCriterion)) #TODO not very robust
        while(not flag['stop']):
            # TODO Mostly outdated. Must be thoroughly checked.
            newHierarchy,splittingParameter = self.optimalHierarchy()

            self.splitTolerance(splittingParameter)

            self.updateHierarchy(newHierarchy)

            # synchronization point needed to launch new tasks if convergence is false
            # put the synchronization point as in the end as possible
            # TODO: remove from here the synchronization point to more hidden places
            flag = self.stoppingCriterionFlag(self.iterationCounter)

            # TODO Print selection is mostly guesswork here (very fragile)
            print("Iteration ",self.iterationCounter, "\tTolerance - ", ['%.3e' % tol for tol in self.tolerances(self.tolerancesForHierarchy) ],
                    "\tError - ", ['%.3e' % err for err in get_value_from_remote(self.errorEstimation(self.errorsForStoppingCriterion))],
                "\tHierarchy - ", self.hierarchy())

            if flag['updateTolerance']:
                self.updateTolerance()
            if flag['updateIndexSpace']:
                self.updateHierarchySpace()

            self.iterationCounter += 1

            #### DATA DUMP ##########
            if self.isDataDumped is True:
                pathObject = pl.Path(self.outputFolderPath)
                pathObject.mkdir(parents=True,exist_ok=True)
                filename = self.outputFolderPath+"/iteration_"+str(self.iterationCounter)+".pickle"
                output_file = open(filename,"wb")

                output_dict = {}
                hier = self.hierarchy()
                output_dict['predictedHierarchy']=newHierarchy
                output_dict['hierarchy']=hier
                if(len(self.predictorsForHierarchy)!=0):
                    qoip = self.qoiPrediction()
                    costp = self.costPrediction()
                    output_dict['biasParameters']=get_value_from_remote(qoip[0].parameters)
                    output_dict['varParameters']=get_value_from_remote(qoip[1].parameters)
                    output_dict['costParameters']=get_value_from_remote(costp.parameters)

                output_dict['indexwiseBias']=get_value_from_remote(self.indexEstimation(0,[1,True,False]))
                errs = get_value_from_remote(self.indexEstimation(0,[1,True,True]))
                levels,samples = splitOneListIntoTwo(hier)
                output_dict['indexwiseVar']=[errs[i]*samples[i] for i in range(len(errs))]
                output_dict['indexwiseCost']=get_value_from_remote(self.indexCostEstimation([1,True,False]))

                hier = newHierarchy
                levels,samples = splitOneListIntoTwo(hier)
                costs = get_value_from_remote(self.indexCostEstimation([1,True,False]))
                total_times = [sample*cost for sample,cost in zip(samples,costs)]
                output_dict['totalTime']=sum(total_times)
                pickle.dump(output_dict,output_file)

        # TODO - Debug statement. Remove for PYCOMPSS tests
        print("Estimation - ",get_value_from_remote(self.estimation()),
            "\nFinal Error - ",get_value_from_remote(self.errorEstimation(self.errorsForStoppingCriterion)),
            "\nLevelwise Mean Times - ",get_value_from_remote(self.indexCostEstimation([1,True,False])))


    ####################################################################################################
    ###################################### ASYNCHRONOUS FRAMEWORK ######################################
    ####################################################################################################


    def asynchronousFinalizeIteration(self):
        continue_iterating = True
        for batch in range (self.monteCarloSampler.numberBatches):
            if (self.monteCarloSampler.batchesLaunched[batch] is True and self.monteCarloSampler.batchesExecutionFinished[batch] is True and self.monteCarloSampler.batchesAnalysisFinished[batch] is True and self.monteCarloSampler.batchesConvergenceFinished[batch] is not True and continue_iterating is True):
                continue_iterating = False
                self.monteCarloSampler.asynchronousFinalize(batch)
                flag = self.stoppingCriterionFlag(self.iterationCounter)
                self.monteCarloSampler.batchesConvergenceFinished[batch] = True
                break
        # screen iteration informations
        print("Iteration ",self.iterationCounter, "\tTolerance - ", ['%.3e' % tol for tol in self.tolerances(self.tolerancesForHierarchy) ],
                "\tError - ", ['%.3e' % err for err in get_value_from_remote(self.errorEstimation(self.errorsForStoppingCriterion))],
            "\tHierarchy - ", self.hierarchy())
        # update tolerance and hierarchy space if required
        if flag['updateTolerance']:
            self.updateTolerance()
        if flag['updateIndexSpace']:
            self.updateHierarchySpace()
        # update iteration counter
        self.iterationCounter += 1
        return flag

    def runAsynchronousXMC(self):
        """
        Run specified algorithm with asynchronous framework.
        """

        self.checkInitialisation()
        # set maximum number of iteration variable
        if (self.stoppingCriterion.tolerances([self.positionMaxNumberIterationsCriterion])[0] is not None):
            if (type(self.stoppingCriterion.tolerances([self.positionMaxNumberIterationsCriterion])[0]) is int):
                if (self.stoppingCriterion.tolerances([self.positionMaxNumberIterationsCriterion])[0] != self.stoppingCriterion.tolerances([-1])[0]):
                    raise Exception ("Set the correct positionMaxNumberIterationsCriterion in XMC algorithm inpuct dictionary. It should be the last entry of monoCriteriaInpuctDict.")
                self.monteCarloSampler.maximumNumberIterations = self.stoppingCriterion.tolerances([self.positionMaxNumberIterationsCriterion])[0]
            else:
                raise Exception ("positionMaxNumberIterationsCriterion not set in XMC algorithm settings. Please set it in order to run the asynchronous framework.")
        else:
            raise Exception ("positionMaxNumberIterationsCriterion not set in XMC algorithm settings. Please set it in order to run the asynchronous framework.")
        # Iteration loop will start here
        flag = self.stoppingCriterion.flagStructure()
        self.iterationCounter = 0
        while(not flag['stop']):
            # estimate splitting parameter
            newHierarchy,splittingParameter = self.optimalHierarchy()
            self.splitTolerance(splittingParameter)
            # launch asynchronous monteCarloSampler update method
            self.monteCarloSampler.asynchronousUpdate(newHierarchy)
            # finalize phase of the iteration and return flag
            flag = self.asynchronousFinalizeIteration()

        # screen results
        print("Estimation - ",get_value_from_remote(self.estimation()),
            "\nFinal Error - ",get_value_from_remote(self.errorEstimation(self.errorsForStoppingCriterion)),
            "\nLevelwise Mean Times - ",get_value_from_remote(self.indexCostEstimation([1,True,False])))
