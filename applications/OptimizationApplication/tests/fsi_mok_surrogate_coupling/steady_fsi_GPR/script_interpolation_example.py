#!/usr/bin/env python
###############################################################################
# General example of using the smarty.interpolation module.
# Samples from an analytical test function and validates the surrogate
# on a generated test data set and via KFold cross-validation.
#
# In addition, the use of the doe, snapshot and the Log module are shown.
# The switch 'doComplex' switch the target values to be complex valued. This
# does not need any changes with the model, only for target value generation.
#
###############################################################################
# Author          : M. Stradtner
# Date of creation: 17/01/2019
###############################################################################
import os
import numpy as np

from smarty import interpolation, Log, snapshot
from smarty.sampling import doe

from smarty.modelSelection import ModelSelection

###############################################################################
# --- OUTPUT TMP FOLDER ---
###############################################################################
# Please set the path for the tmp output to the desired folder or change the
# output directory in the USER INPUT.
# By default, the application creates a tmp folder in the SMARTY APP folder.
###############################################################################
from applications import SMARTY_DATA_PATH, SMARTY_TMP
OUTPUT_PATH = SMARTY_TMP

###############################################################################
# --- USER INPUT ---
###############################################################################
userInput = {
    "doeMethod"                 : doe.Halton,
    "doeOptions"                : {"designSpace"        : [(-3., 3), (-3., 3.)],
                                   "numberOfSamples"    : 20,
                                   "parameterNames"     : ["x", "y"],
                                   },
    "interpolationMethod"       : interpolation.Gaussian,  # interpolation.Kriging, interpolation.TPS
    "interpolationOptions"      : {"augmentation"       : 1,
                                   "regularization"     : True,
                                   "lowFidelityModel"   : None,
                                   "tuneModus"          : 'ML',  # None, 'CV', 'RandomSearch'
                                   "tuneOptions"        : None,
                                   "optimizationOptions": None,
                                   },
    "errorMetrics"              : ["root mean squared error", "maximum error", "r2 score"],  # defaults
    "outputDir"                 : "/home/virtro/Desktop/Tutorials-master/multiphysics/steady_fsi",
    "nFoldsCV"                  : 5,
    "doComplex"                 : False,
    "nRepeats"                  : 3
}


def testfun(x):
    return (3 * (1 - x[0]) ** 2 * np.exp(-1.0 * x[0] ** 2 - (x[1] + 1) ** 2)
            - 10 * (x[0] / 5 - x[0] ** 3 - x[1] ** 5) * np.exp(-x[0] ** 2 - x[1] ** 2)
            - 1 / 3 * np.exp(-(x[0] + 1) ** 2 - x[1] ** 2))


##############################################################################
# --- MAIN FUNCTION ---
##############################################################################
def main():
    # Create output folder
    if not os.path.exists(userInput["outputDir"]):
        Log("Create output directory: {0}".format(userInput["outputDir"]))
        os.makedirs(userInput["outputDir"])
    else:
        Log("Output directory: {0}".format(userInput["outputDir"]))

    # Get a data sets (training and test set).
    Log("Get samples ...")
    trainingSamples = getSamples("TrainingSet", userInput["doComplex"],
                                 userInput["doeMethod"], **userInput["doeOptions"])
    newSamples = getSamples("NewSamples", userInput["doComplex"], doe.FullFactorial,
                            **{"parameters": [(-3, 3, 50), (-3, 3, 50)], "parameterNames": ["x", "y"]})
    testSamples = getSamples("TestSet", userInput["doComplex"], doe.FullFactorial,
                             **{"parameters": [(-3, 3, 10), (-3, 3, 10)], "parameterNames": ["x", "y"]})
    Log("Snapshots created.")

    # Build a surrogate model.
    Log("Build a surrogate model ...")
    surrogate = buildSurrogate(trainingSamples,
                               interpolationMethod=userInput["interpolationMethod"],
                               interpolationOptions=userInput["interpolationOptions"])
    Log("Surrogate model built.")

    outputFilename = os.path.join(userInput["outputDir"], "predictions.nc")
    Log("Predict samples and write to file {0}".format(outputFilename))
    predictions = getPredictions(newSamples, surrogate)
    writeOutput(predictions, outputFilename, asciiTecplotFormat=True)

    # validate your surrogate model
    # on given test data
    validate(testSamples, surrogate, errorMetrics=userInput["errorMetrics"])

    # via cross-validation
    kwargsCV = userInput["interpolationOptions"].copy()
    kwargsCV["interpolationMethod"] = userInput["interpolationMethod"]
    crossValidation(trainingSamples,
                    surrogateParameter=kwargsCV,
                    errorMetrics=userInput["errorMetrics"],
                    kFolds=userInput["nFoldsCV"],
                    nRepeats=userInput["nRepeats"])

    Log("Main()::Done.")


###############################################################################
# --- SUPPORT FUNCTIONS ---
###############################################################################
def getSamples(case, doComplex, doeMethod, **doeOptions):
    """
    Computes a sampling plan (Design of Experiment), evaluates the cost
    function for the given sample locations and creates a SMARTy snapshot for
    data handling.
    """
    sampling = doeMethod(**doeOptions)
    controlPoints = sampling.GetSampling()

    if doComplex:
        numtype = np.complex
    else:
        numtype = float
    values = np.zeros((controlPoints.shape[0], 1), dtype=numtype)
    for i in range(controlPoints.shape[0]):
        values[i] = testfun(controlPoints[i])
        if doComplex:
            values[i] += 1j * testfun(controlPoints[-i])

    parameters = {"dataset": case}
    variables = {"x": controlPoints[:, 0],
                 "y": controlPoints[:, 1],
                 "value": values}

    return snapshot.Snapshot(variables=variables, attributes=parameters)


def buildSurrogate(snaps, interpolationMethod, interpolationOptions=None):
    """
    Builds a surrogate model (TPS, Gaussian, etc.)
    """
    interpolationOptions = {} if interpolationOptions is None else interpolationOptions
    options = interpolationOptions.copy()
    tuneModus = options.pop("tuneModus", None)
    tuneOptions = options.pop("tuneOptions", None)
    optimizationOptions = options.pop("optimizationOptions", None)

    controlPoints = snaps.GetValues(["x", "y"], asMatrix=True)
    values = snaps.GetValues(["value"], asMatrix=True)

    return interpolationMethod(**options).Fit(controlPoints, values, tuneModus=tuneModus,
                                              tuneOptions=tuneOptions,
                                              optimizationOptions=optimizationOptions)


def getPredictions(snaps, surrogate):
    """
    Predict new samples at sample sites specified by the given snapshot.
    """
    controlPoints = snaps.GetValues(["x", "y"], asMatrix=True)
    values = surrogate.GetPrediction(controlPoints)

    snaps._Variables.update({"values": values})

    return snaps


def validate(snaps, surrogate, errorMetrics):
    """
    Validates the surrogate model
    """
    controlPoints = snaps.GetValues(["x", "y"], asMatrix=True)
    testValues = snaps.GetValues(["value"], asMatrix=True)

    errors = surrogate.Validate(controlPoints, testValues, errorMetrics=errorMetrics)

    for k in errors:
        Log("{0}: {1}".format(k, errors[k]))


def crossValidation(snaps, surrogateParameter, errorMetrics, kFolds=5, nRepeats=5):
    """
    Run a KFold Cross-Validation.
    """
    controlPoints = snaps.GetValues(["x", "y"], asMatrix=True)
    values = snaps.GetValues(["value"], asMatrix=True)

    ms = ModelSelection(controlPoints, values, surrogateParameter, cvKFolds=kFolds,
                        errorMetrics=errorMetrics, nRepeats=nRepeats)

    ms.RunModelSelection()
    errors = ms.GetModelSelectionResults()

    for k in sorted(errors.keys()):
        Log("\n{0}: {1}".format(k, errors[k]))


def writeOutput(snaps, outputFilename, asciiTecplotFormat=False):
    """
    Writes the Snapshot to file (NetCDF format). Into ASCII-Tecplot format if
    option enabled.
    """
    snaps.WriteToFile(outputFilename)
    if asciiTecplotFormat:
        outputFilename = os.path.splitext(outputFilename)[0] + ".dat"
        variables = ["x", "y", "value"]
        snaps.WriteToTECPLOTFile(outputFilename, variables, fmt='+10.6f')


###############################################################################
if __name__ == '__main__':
    main()
