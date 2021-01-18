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
from exaqute.ExaquteTask import *

from pycompss.api.mpi import mpi 
from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.api import compss_barrier
from pycompss.api.api import compss_delete_object
from pycompss.api.api import compss_delete_file

from pycompss.api.parameter import *

from pycompss.api.implement import implement

from pycompss.api.constraint import *



#class ExaquteTask(object):
#
#    def __init__(self, *args, **kwargs):
#        global scheduler
#        scheduler = "Current scheduler is PyCOMPSs"
#        self.task_instance = task(*args, **kwargs)
#
#    def __call__(self, f):
#        return self.task_instance.__call__(f)

ExaquteTask = task


def barrier():  # Wait
    compss_barrier()


def get_value_from_remote(obj):  # Gather
    obj = compss_wait_on(obj)
    return obj


def delete_object(*objs):  # Release
    for obj in objs:
        compss_delete_object(obj)


def delete_file(file_path):
    compss_delete_file(file_path)


def compute(obj):  # Submit task
    return obj
