from exaqute.ExaquteTask import *

from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.api import compss_barrier
from pycompss.api.api import compss_delete_object
from pycompss.api.api import compss_delete_file

from pycompss.api.parameter import *

from pycompss.api.implement import implement

from pycompss.api.constraint import *


class ExaquteTask(object):

    def __init__(self, *args, **kwargs):
        global scheduler
        scheduler = "Current scheduler is PyCOMPSs"
        self.task_instance = task(*args, **kwargs)

    def __call__(self, f):
        return self.task_instance.__call__(f)


def barrier():  # Wait
    compss_barrier()


def get_value_from_remote(obj):  # Gather
    obj = compss_wait_on(obj)
    return obj


def delete_object(obj):  # Release
    compss_delete_object(obj)


def delete_file(file_path):
    compss_delete_file(file_path)


def compute(obj):  # Submit task
    return obj
