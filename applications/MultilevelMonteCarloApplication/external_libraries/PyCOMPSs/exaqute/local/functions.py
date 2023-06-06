from exaqute.common import ExaquteException

from .internals import (
    ValueWrapper,
    _check_accessed,
    _check_init,
    _delete_accessed,
    _init,
    _submit_point,
)


def init():
    _init()


def get_value_from_remote(obj):
    _check_init()
    _submit_point()
    t = type(obj)
    if t is list:
        return [get_value_from_remote(o) for o in obj]
    if t is ValueWrapper:
        if not obj.keep:
            raise ExaquteException(
                "get_value_from_remote called on not keeped object, object created at {}".format(
                    obj.traceback
                )
            )
        return obj.unwrap_value()
    else:
        if isinstance(obj, (int, bool, float, str)):
            return obj
        else:
            if not _check_accessed(obj):
                print(
                    "WARN: get_value_from_remote called on non-task value, got {!r}".format(
                        obj
                    )
                )
            else:
                _delete_accessed(obj)
            return obj


def barrier():
    _check_init()
    _submit_point()


def delete_object(*objs):
    for obj in objs:
        _delete_object(obj)


def _delete_object(obj):
    _check_init()
    t = type(obj)
    if t is list:
        for o in obj:
            _delete_object(o)
    elif t is ValueWrapper:
        if not obj.keep:
            raise ExaquteException(
                "Deleting non-keeped object, object created at {}".format(obj.traceback)
            )
        if obj.deleted:
            raise ExaquteException(
                "Deleting already deleted object, object created at {}".format(
                    obj.traceback
                )
            )
        obj.deleted = True
    else:
        if not isinstance(obj, (bool, int, float, str)):
            if not _check_accessed(obj):
                raise ExaquteException(
                    "delete_object called on non-task value, got {!r}".format(obj)
                    + " type "
                    + str(t)
                )
            else:
                _delete_accessed(obj)
