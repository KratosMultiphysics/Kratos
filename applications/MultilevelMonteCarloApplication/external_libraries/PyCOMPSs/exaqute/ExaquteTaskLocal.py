from exaqute.ExaquteTask import *


class ExaquteTask(object):

    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, f):
        def g(*args, **kwargs):
            if "scheduling_constraints" in kwargs:
                del kwargs["scheduling_constraints"]
            return f(*args, **kwargs)

        return g


def barrier():  # Wait
    pass


def get_value_from_remote(obj):  # Gather
    return obj


def delete_object(obj):  # Release
    del obj


def delete_file(file_path):
    import os
    os.remove(file_path)


def compute(obj):  # Submit task
    return obj


class implement(object):

    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, f):
        return f


class constraint(object):

    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, f):
        return f
