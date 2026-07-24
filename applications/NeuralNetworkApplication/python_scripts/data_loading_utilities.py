import numpy as np
import h5py
import ast
import json
import pickle
from KratosMultiphysics.NeuralNetworkApplication.input_dataclasses import (
    ListDataWithLookback,
    ListNeuralNetworkData,
)


def ImportDataFromFile(external_file, category, lookback=False):
    "Import the data from a file generated for neural network training."
    # Data loading for h5
    if external_file.endswith(".h5"):
        data = ImportH5(external_file, lookback)
    # Data loading for dat
    elif external_file.endswith(".dat") or external_file.endswith(".csv"):
        data = ImportAscii(external_file, lookback)
    # Data loading for npy
    elif external_file.endswith(".npy"):
        data = ImportNpy(external_file, lookback)
    # Data loading for pkl
    elif external_file.endswith(".pkl"):
        data = ImportPkl(external_file)
    # Exception for non-supported formats
    else:
        raise Exception(
            category
            + " data format not supported. Supported formats are .dat, .npy, .pkl and .h5"
        )
    return data


def ImportH5(external_file, lookback):
    "Import the data of an h5 file generated for neural network training. Check the expected format."
    print("Warning: HDF5 functionalities are not fully implemented.")
    if lookback:
        raw = ListDataWithLookback()
    else:
        raw = ListNeuralNetworkData()
    with h5py.File(external_file, "r") as f:
        for time in f.keys():
            time_array = []
            for model_part in f[time].keys():
                for variables in f[time][model_part]["NodalSolutionStepData"].keys():
                    value = f[time][model_part]["NodalSolutionStepData"][variables][:]
                    time_array.append(value[:])
            raw.AddToList(np.squeeze(time_array))
    return raw


def ImportAscii(external_file, lookback):
    "Import the data of an ascii file generated for neural network training."
    raw = []
    with open(external_file, "r") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            raw.append(list(map(str, line.split())))
    if lookback:
        raw_splits = ListDataWithLookback()
    else:
        raw_splits = ListNeuralNetworkData()
    for line in raw:
        raw_splits_line = []
        for data in line:
            if not isfloat(data.split(" ")[0]):
                array = "".join(data[4:-1])
                array_floats = list(
                    map(
                        float,
                        array.split(
                            ",",
                        ),
                    )
                )
                raw_splits_line.extend(array_floats)
            else:
                array_floats = float(data)
                raw_splits_line.append(array_floats)
        raw_splits.AddToList(np.squeeze(raw_splits_line))
    return raw_splits


def ImportNpy(external_file, lookback):
    "Import the data from a npy file."
    if lookback:
        raw = ListDataWithLookback()
    else:
        raw = ListNeuralNetworkData()
    raw.ImportFromArray(np.squeeze(np.load(external_file)))
    return raw


def ImportPkl(external_file):
    "Import the data from a pkl file that contains a ListNeurakNetworkData or ListDataWithLookback object."
    with open(external_file, "rb") as file:
        raw = pickle.load(file)
    return raw


def ImportDictionaryFromText(external_file):
    "Import the data from a text file to a dictionary"

    with open(external_file, "r") as file:
        data = file.read()

    return ast.literal_eval(data)


def UpdateDictionaryJson(external_file, dump_data):
    with open(external_file, "w") as f:
        f.seek(0)
        json.dump(dump_data, f)
        f.truncate()


def KratosVectorToList(vector):
    vectortolist = []
    for i in vector:
        vectortolist.append(i)
    return vectortolist


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False
