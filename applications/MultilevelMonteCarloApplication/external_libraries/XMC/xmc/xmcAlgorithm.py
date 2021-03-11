import pickle
import pathlib as pl

# XMC imports
from xmc.tools import dynamicImport, splitOneListIntoTwo
from xmc.methodDefs_xmcAlgorithm import checkInitialisation, updateTolerance
from xmc.distributedEnvironmentFramework import get_value_from_remote


class XMCAlgorithm:
    """
    This top-level class handles the overall algorithm: initialisation as well as everything
    related to error, tolerance and iterations. It also possesses the necessary methods and
    attributes to create new types of algorithms. However, the export of results is to be
    handled outside.
    This documentation is not yet complete.

    * Methods
    - runXMC: run an algorithm with generic structure. It is the method to call to run the
      algorithm.
    - runAsynchronousXMC: run an algorithm with generic structure exploiting the asynchronous
      framework. It is the method to call to run the asynchronous algorithm.

    * Attributes
    - estimatorsForHierarchy: List[List]. This is a list of instructions for the indexwise
    estimations to be sent to the hierarchy optimiser. estimatorsForHierarchy[0] contains
    instructions for the first estimation to be sent (estimatorsForHierarchy[1] the second,
    etc.). estimatorsForHierarchy[0][0] is the index of the estimator concerned, as ordered in
    MonteCarloIndex.qoiEstimator; estimatorsForHierarchy[0][1] is the list of arguments to be
    passed to the value method of this estimator. These values are used in
    XMCAlgorithm.optimalHierarchy and eventually result in a call of the form
    StatisticalEstimator.value(estimatorsForHierarchy[0], *estimatorsForHierarchy[1]).

    """

    def __init__(self, **keywordArgs):
        # TODO - think of a better way to do data dumps and output
        self.isDataDumped = keywordArgs.get("isDataDumped", False)
        self.outputFolderPath = keywordArgs.get("outputFolderPath", None)

        # Attributes
        self.monteCarloSampler = keywordArgs.get("monteCarloSampler")
        self.hierarchyOptimiser = keywordArgs.get("hierarchyOptimiser")
        self.stoppingCriterion = keywordArgs.get("stoppingCriterion")
        self.errorsForStoppingCriterion = keywordArgs.get("errorForStoppingCriterion", None)
        self.predictorsForHierarchy = keywordArgs.get("predictorsForHierarchy", None)
        self.estimatorsForHierarchy = keywordArgs.get("estimatorsForHierarchy", None)
        self.tolerancesForHierarchy = keywordArgs.get("tolerancesForHierarchy", None)
        self.errorParametersForHierarchy = keywordArgs.get("errorParametersForHierarchy", None)
        self.costPredictorForHierarchy = keywordArgs.get("costPredictorForHierarchy", None)
        self.costEstimatorForHierarchy = keywordArgs.get("costEstimatorForHierarchy", None)
        self.positionMaxNumberIterationsCriterion = keywordArgs.get(
            "positionMaxNumberIterationsCriterion", None
        )
        self.iterationCounter = 0

        # Methods
        self.checkInitialisation = dynamicImport(
            keywordArgs.get("checkInitialisation", "xmc.tools.doNothing")
        )

    def updateTolerance(self, criteriaToUpdate=None):
        self.stoppingCriterion.updateTolerance(criteriaToUpdate)

    def indexEstimation(self, coordinate, valueMethodArgs):
        return self.monteCarloSampler.indexEstimation(coordinate, valueMethodArgs)

    def indexCostEstimation(self, valueMethodArgs):
        return self.monteCarloSampler.indexCostEstimation(valueMethodArgs)

    def predictor(self, choice=None):
        if choice is None:
            choice = slice(0, len(self.monteCarloSampler.qoiPredictor))
        return self.monteCarloSampler.qoiPredictor[choice]

    def costPredictor(self):
        return self.monteCarloSampler.costPredictor

    def tolerances(self, *args):
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
        input_dict = self.hierarchyOptimiser.inputDictionaryTemplate()

        # No optimisation at first iteration
        if self.iterationCounter < 1:
            newHierarchy = self.hierarchyOptimiser.defaultHierarchy
            splittingParameter = input_dict.get("splittingParameter", None)
            return newHierarchy, splittingParameter

        # Else, assemble data for hierarchy optimiser

        # Random variables of interest
        # Indexwise estimations
        input_dict["estimations"] = [
            self.indexEstimation(c[0], c[1]) for c in self.estimatorsForHierarchy
        ]
        # Predictors
        if self.predictorsForHierarchy:
            input_dict["models"] = []
            input_dict["parametersForModel"] = []
            for coord in self.predictorsForHierarchy:
                input_dict["models"].append(self.predictor(coord)._valueForParameters)
                # TODO This should get self.predictor(coord).oldParameters
                # and default to self.predictor(coord).parameters if they are None
                input_dict["parametersForModel"].append(self.predictor(coord).parameters)

        # Sample cost
        # Indexwise estimation
        if self.costEstimatorForHierarchy is not None:
            input_dict["costEstimations"] = self.indexCostEstimation(
                self.costEstimatorForHierarchy
            )

        # Predictor
        if self.costPredictor() is not None:
            input_dict["costModel"] = self.costPredictor()._valueForParameters
            # TODO This should get self.costPredictor().oldParameters
            # and default to self.costPredictor().parameters if they are None
            input_dict["costParameters"] = self.costPredictor().parameters

        # Error parameters
        # TODO - Triple dereference below!! Add method to get errorEstimator parameters
        # or errorEstimator objects themselves from monteCarloSampler
        if self.errorParametersForHierarchy is not None:
            input_dict["errorParameters"] = [
                self.monteCarloSampler.errorEstimators[c].parameters
                for c in self.errorParametersForHierarchy
            ]

        # Miscellaneous parameters
        input_dict["newSampleNumber"] = 25  # TODO configurable, not hard-coded
        input_dict["oldHierarchy"] = self.hierarchy()
        input_dict["tolerances"] = self.tolerances(self.tolerancesForHierarchy)
        input_dict["defaultHierarchy"] = self.hierarchyOptimiser.defaultHierarchy

        # Synchronisation
        input_dict = get_value_from_remote(input_dict)

        # Compute new hierarchy
        newHierarchy = self.hierarchyOptimiser.optimalHierarchy(input_dict)
        splittingParameter = input_dict.get("splittingParameter", None)
        return newHierarchy, splittingParameter

    def updateHierarchy(self, newHierarchy):
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

    def updateHierarchySpace(self, *args):
        """
        Method that interfaces with the HierarchyOptimiser class to compute
        the hierarchy space in which to search for the optimal hierarchy
        """
        self.hierarchyOptimiser.updateHierarchySpace(args)

    def stoppingCriterionFlag(self, currentCost=None):
        """
        Call stoppingCriterion.flag with the proper arguments and return its output
        (a.k.a flag).
        Input argument: currentCost is an indication of the cost the algorithm has entailed
        so far; we usually use the number of iterations.
        Output argument: criterion flag structure as define in the MultiCriterion class.
        """
        # Get errors required for stopping criterion
        errors = self.errorEstimation(self.errorsForStoppingCriterion)

        # Set up dictionary required for stoppingCriterion.flag
        input_dictionary = {}
        for i in range(len(errors)):
            input_dictionary["error" + str(i)] = get_value_from_remote(errors[i])
        input_dictionary["hierarchy"] = self.hierarchy()
        input_dictionary["algorithmCost"] = currentCost

        # Compute flag from dictionary and return
        flag = self.stoppingCriterion.flag(input_dictionary)
        return flag

    def runXMC(self):
        """
        Run an algorithm with generic structure.
        """

        self.checkInitialisation(self)

        # Iteration Loop will start here
        flag = self.stoppingCriterion.flagStructure()
        self.iterationCounter = 0
        # print("Beginning Iterations for tolerance ", self.tolerances(self.errorsForStoppingCriterion)) #TODO not very robust
        while not flag["stop"]:
            # TODO Mostly outdated. Must be thoroughly checked.
            newHierarchy, splittingParameter = self.optimalHierarchy()

            self.splitTolerance(splittingParameter)

            self.updateHierarchy(newHierarchy)

            # synchronization point needed to launch new tasks if convergence is false
            # put the synchronization point as in the end as possible
            # TODO: remove from here the synchronization point to more hidden places
            flag = self.stoppingCriterionFlag(self.iterationCounter)
            flag = get_value_from_remote(flag)

            # TODO Display selection is mostly guesswork here (very fragile)
            errors = get_value_from_remote(
                self.errorEstimation(self.errorsForStoppingCriterion)
            )
            dErrors = " ".join(["{err:.3e}".format(err=float(error)) for error in errors])
            dHierarchy = " ".join([str(i[1]) for i in self.hierarchy()])
            dTol = " ".join(
                [
                    "{t:.3e}".format(t=tol)
                    for tol in self.tolerances(self.tolerancesForHierarchy)
                ]
            )
            print(
                f"Iteration — {self.iterationCounter}",
                f"Tolerances — {dTol}",
                f"Errors — {dErrors}",
                f"Hierarchy — {dHierarchy}",
                sep="\t",
            )

            if flag["updateTolerance"]:
                self.updateTolerance()
            if flag["updateIndexSpace"]:
                self.updateHierarchySpace()

            self.iterationCounter += 1

            #### DATA DUMP ##########
            if self.isDataDumped is True:
                pathObject = pl.Path(self.outputFolderPath)
                pathObject.mkdir(parents=True, exist_ok=True)
                filename = (
                    self.outputFolderPath
                    + "/iteration_"
                    + str(self.iterationCounter)
                    + ".pickle"
                )
                output_file = open(filename, "wb")

                output_dict = {}
                hier = self.hierarchy()
                output_dict["predictedHierarchy"] = newHierarchy
                output_dict["hierarchy"] = hier
                if len(self.predictorsForHierarchy) != 0:
                    qoip = self.predictor()
                    costp = self.costPredictor()
                    output_dict["biasParameters"] = qoip[0].parameters
                    output_dict["varParameters"] = qoip[1].parameters
                    output_dict["costParameters"] = costp.parameters

                output_dict["indexwiseBias"] = self.indexEstimation(0, [1, True, False])
                errs = self.indexEstimation(0, [1, True, True])
                levels, samples = splitOneListIntoTwo(hier)
                output_dict["indexwiseVar"] = [errs[i] * samples[i] for i in range(len(errs))]
                output_dict["indexwiseCost"] = self.indexCostEstimation([1, True, False])

                hier = newHierarchy
                levels, samples = splitOneListIntoTwo(hier)
                costs = self.indexCostEstimation([1, True, False])
                total_times = [sample * cost for sample, cost in zip(samples, costs)]
                output_dict["totalTime"] = sum(total_times)
                pickle.dump(output_dict, output_file)

        # TODO - Debug statement. Remove for PYCOMPSS tests
        displayEstimation = get_value_from_remote(self.estimation())
        displayEstimation = " ".join(["{e:.3e}".format(e=est) for est in displayEstimation])
        displayError = get_value_from_remote(
            self.errorEstimation(self.errorsForStoppingCriterion)
        )
        displayError = " ".join(["{e:.3e}".format(e=error) for error in displayError])
        displayCost = get_value_from_remote(self.indexCostEstimation([1, True, False]))
        displayCost = " ".join(["{c:.3e}".format(c=cost) for cost in displayCost])
        print(
            f"Estimations — {displayEstimation}",
            f"Final errors — {displayError}",
            f"Levelwise mean costs — {displayCost}",
            sep="\n",
        )

    ####################################################################################################
    ###################################### ASYNCHRONOUS FRAMEWORK ######################################
    ####################################################################################################

    def asynchronousFinalizeIteration(self):
        """
        Method finalizing an iteration of the asynchornous framework. It synchronizes and calls all relevant methods of one single batch, the first available, before estimating convergence.
        """

        continue_iterating = True
        for batch in range(self.monteCarloSampler.numberBatches):
            if (
                self.monteCarloSampler.batchesLaunched[batch] is True
                and self.monteCarloSampler.batchesExecutionFinished[batch] is True
                and self.monteCarloSampler.batchesAnalysisFinished[batch] is True
                and self.monteCarloSampler.batchesConvergenceFinished[batch] is not True
                and continue_iterating is True
            ):
                continue_iterating = False
                self.monteCarloSampler.asynchronousFinalize(batch)
                flag = self.stoppingCriterionFlag(self.iterationCounter)
                self.monteCarloSampler.batchesConvergenceFinished[batch] = True
                break
        # screen iteration informations
        errors = get_value_from_remote(self.errorEstimation(self.errorsForStoppingCriterion))
        print(
            "Iteration ",
            self.iterationCounter,
            "\tTolerance - ",
            ["%.3e" % tol for tol in self.tolerances(self.tolerancesForHierarchy)],
            "\tError - ",
            ["%.3e" % err for err in errors],
            "\tHierarchy - ",
            self.hierarchy(),
        )
        # update tolerance and hierarchy space if required
        if flag["updateTolerance"]:
            self.updateTolerance()
        if flag["updateIndexSpace"]:
            self.updateHierarchySpace()
        # update iteration counter
        self.iterationCounter += 1
        return flag

    def runAsynchronousXMC(self):
        """
        Run algorithm with asynchronous framework.
        """
        self.checkInitialisation(self)
        # set maximum number of iteration variable
        self.monteCarloSampler.maximumNumberIterations = self.stoppingCriterion.tolerances(
            [self.positionMaxNumberIterationsCriterion]
        )[0]
        # Iteration loop will start here
        flag = self.stoppingCriterion.flagStructure()
        self.iterationCounter = 0
        while not flag["stop"]:
            # estimate splitting parameter
            newHierarchy, splittingParameter = self.optimalHierarchy()
            self.splitTolerance(splittingParameter)
            # launch asynchronous monteCarloSampler update method
            self.monteCarloSampler.asynchronousUpdate(newHierarchy)
            # finalize phase of the iteration and return flag
            flag = self.asynchronousFinalizeIteration()

        # screen results
        displayEstimation = self.estimation()
        displayError = self.errorEstimation(self.errorsForStoppingCriterion)
        displayCost = self.indexCostEstimation([1, True, False])
        print(
            f"Estimations — {displayEstimation}",
            f"Final errors — {displayError}",
            f"Levelwise mean costs — {displayCost}",
            sep="\n",
        )
