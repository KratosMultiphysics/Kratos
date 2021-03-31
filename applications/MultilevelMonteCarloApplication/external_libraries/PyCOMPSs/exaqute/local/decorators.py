import os
from functools import wraps

from exaqute.common import ExaquteException

from .consts import (
    COLLECTION_IN,
    COLLECTION_INOUT,
    COLLECTION_OUT,
    FILE_IN,
    FILE_INOUT,
    FILE_OUT,
    IN,
    INOUT,
)
from .internals import ValueWrapper, _check_init, _is_nested, _obj_to_value


class Task(object):
    """
    Dummy task class (decorator style)
    """

    def __init__(self, *args, **kwargs):
        self.args = args
        self.returns = 1
        self.priority = False
        self.keep = False
        if "returns" in kwargs:
            self.returns = kwargs.pop("returns")
        if "priority" in kwargs:
            self.priority = kwargs.pop("priority")
        if "keep" in kwargs:
            self.keep = kwargs.pop("keep")
        for k, v in kwargs.items():
            if v not in (
                IN,
                INOUT,
                FILE_IN,
                FILE_OUT,
                FILE_INOUT,
                COLLECTION_IN,
                COLLECTION_OUT,
                COLLECTION_INOUT,
            ):
                if not isinstance(v, dict):
                    raise ExaquteException("Argument '{}' has invalid value".format(k))

    def __call__(self, f):
        @wraps(f)
        def wrapped_task(*args, **kwargs):
            _check_init()
            if _is_nested():
                return f(*args, **kwargs)
            if "returns" in kwargs:
                returns = kwargs.pop("returns")
            else:
                returns = self.returns
            if "keep" in kwargs:
                keep = kwargs.pop("keep")
            else:
                keep = self.keep
            new_args = [_obj_to_value(obj) for obj in args]
            for k, v in kwargs.items():
                kwargs[k] = _obj_to_value(v)
            result = f(*new_args, **kwargs)
            if returns != 1:
                if not hasattr(result, "__len__") or len(result) != returns:
                    raise ExaquteException(
                        "Invalid number of results returned (expected {})".format(
                            returns
                        )
                    )
                return [ValueWrapper(r, keep) for r in result]
            else:
                return ValueWrapper(result, keep)

        return wrapped_task


task = Task


def constraint(computing_units=1):
    if isinstance(computing_units, str):
        try:
            # can be cast to int
            computing_units = int(computing_units)
        except Exception:
            if computing_units.startswith("$"):
                env_var = (
                    computing_units.replace("$", "").replace("{", "").replace("}", "")
                )
                try:
                    computing_units = int(os.environ[env_var])
                except Exception:
                    raise ExaquteException(
                        "Environment var: "
                        + env_var
                        + " not defined or can't be cast to int "
                    )
    assert isinstance(computing_units, int) and computing_units > 0
    return lambda x: x


class Mpi(object):
    def __init__(self, *args, **kwargs):
        self.processes = 1
        if "processes" in kwargs:
            self.processes = kwargs["processes"]
        else:
            raise ExaquteException("Number of processes not specified")
        if "runner" not in kwargs:
            raise ExaquteException("Runner must be specified")
        for k, v in kwargs.items():
            if k not in ("processes", "runner", "flags"):
                if not k.endswith("_layout"):
                    raise ExaquteException("Argument '{}' is not valid".format(k))

    def __call__(self, f):
        def wrapped_f(*args, **kwargs):
            _check_init()
            return f(*args, **kwargs)

        return wrapped_f


mpi = Mpi
