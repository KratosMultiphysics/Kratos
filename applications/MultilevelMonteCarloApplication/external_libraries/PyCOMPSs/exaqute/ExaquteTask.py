#!/usr/bin/python
#
#  Copyright 2002-2019 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
from exaqute.ExaquteParameter import *

class mpi(object):
    def __init__(self, *args, **kwargs):
        raise Exception(
            "Exaqute MPI decorator not implemented in the current scheduler")

    def __call__(self, f):
        raise Exception(
            "Exaqute MPI call code not implemented in the current scheduler")

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
