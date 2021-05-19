import numpy as np
import h5py
import ast
import json


def ImportH5(external_file,category):
    "Import the data of an h5 file generated for neural network training."
    raw = []
    with h5py.File(external_file,'r') as f:
        for time in f.keys():
            time_array = []
            for variables in f[time][category]['NodalSolutionStepData'].keys():
                value = f[time][category]['NodalSolutionStepData'][variables][:]
                # time_array.append(value[:])
                time_array = value[:]
            raw.append(time_array)
    return np.array(raw)

def ImportAscii(external_file):
    "Import the data of an ascii file generated for neural network training."
    raw = []
    with open(external_file,'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            raw.append(list(map(str, line.split())))
    raw_splits =[]
    for line in raw:
        raw_splits_line = []
        for data in line:
            if not isfloat(data.split(' ')[0]):
                array = ''.join(data[4:-1])
                array_floats = list(map(float,array.split(',',)))
                raw_splits_line.append(array_floats)
            else:
                array_floats = float(data)
                raw_splits_line.append(array_floats)
        raw_splits.append(raw_splits_line)
    raw_splits = np.array(raw_splits)
    return raw_splits

def ImportDictionaryFromText(external_file):
    "Import the data from a text file to a dictionary"
    
    with open(external_file, 'r') as file:
        data = file.read()

    return ast.literal_eval(data)

def UpdateDictionaryJson(external_file, dump_data):
    with open(external_file,'w') as f:
                f.seek(0)
                json.dump(dump_data,f)
                f.truncate()
                
def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False