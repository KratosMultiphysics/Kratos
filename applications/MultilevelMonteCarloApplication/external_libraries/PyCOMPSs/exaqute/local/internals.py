import traceback

from exaqute.common import ExaquteException

_exaqute_inited = False
_temp_objects = None
_accessed_objs = []


def _traceback(pop_elems=1):
    stack = traceback.extract_stack(None, 6)
    stack.pop()  # Entry is not usefull, it is actually line above
    return "".join(traceback.format_list(stack))


def _is_nested():
    stack = traceback.extract_stack()
    stack.pop()
    stack.pop()
    st = "".join(traceback.format_list(stack))
    nested = "wrapped_task" in st
    return nested


class ValueWrapper:
    def __init__(self, value, keep):
        self.value = value
        self.keep = keep
        self.deleted = False
        self.traceback = _traceback()
        if not keep:
            _register_temp_object(self)

    def unwrap_value(self):
        if self.deleted:
            raise ExaquteException("Using deleted object")
        if not self.keep and self not in _temp_objects:
            raise ExaquteException(
                "Using temporary object after submit point, object created at {}",
                self.traceback,
            )
        return _obj_to_value(self.value)


def _obj_to_value(obj):
    t = type(obj)
    if t is list:
        return [_obj_to_value(o) for o in obj]
    if t is ValueWrapper:
        return obj.unwrap_value()
    else:
        if isinstance(obj, (int, bool, float, str)):
            return obj
        else:
            if not id(obj) in _accessed_objs:
                _accessed_objs.append(id(obj))
            return obj


def _check_accessed(obj):
    return id(obj) in _accessed_objs


def _delete_accessed(obj):
    _accessed_objs.remove(id(obj))


def _init():
    global _exaqute_inited
    global _temp_objects
    if _exaqute_inited:
        raise ExaquteException("Init called twice")
    _exaqute_inited = True
    _temp_objects = set()


def _reset():
    # For testing purpose
    global _exaqute_inited
    global _temp_objects
    _exaqute_inited = False
    _temp_objects = None


def _register_temp_object(obj):
    _temp_objects.add(obj)


def _check_init():
    if not _exaqute_inited:
        raise ExaquteException("Exaqute call before init")


def _submit_point():
    global _temp_objects
    _temp_objects = set()
