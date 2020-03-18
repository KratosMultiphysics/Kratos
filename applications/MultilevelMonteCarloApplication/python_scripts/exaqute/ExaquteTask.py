from exaqute.ExaquteParameter import *


class ExaquteTask(object):

    def __init__(self, *args, **kwargs):
        raise Exception(
            "Exaqute task decorator not implemented in the current scheduler")

    def __call__(self, f):
        raise Exception(
            "Exaqute task call code not implemented in the current scheduler")


def get_value_from_remote(obj):
    raise Exception("Get value not implemented in the current scheduler")


def barrier():
    raise Exception("Barrier not implemented in the current scheduler")


def delete_object(obj):
    raise Exception("Delete object not implemented in the current scheduler")


def compute(obj):
    raise Exception("Compute not implemented in the current scheduler")

def from_args_to_vector(obj):
    print("FROM ARGS TO VECTOR")
    index = [i for i, x in enumerate(obj) if x == "#"]
    print(index)
    print(list(enumerate(obj)))
    if len(index) == 0:
        return obj
    new_vector = []
    new_vector.append(list(obj[0:index[0]]))
    for i in range(len(index) - 1):
        new_vector.append(list(obj[(index[i] + 1):index[i + 1]]))
    new_vector.append(list(obj[(index[-1] + 1):len(obj)]))
    print("Return:")
    print(new_vector)
    return new_vector


def from_vector_to_args(obj):
    new_vector = []
    new_vector.extend(obj[0])
    for i in range(1, len(obj)):
        new_vector.append("#")
        new_vector.extend(obj[1])
    return new_vector