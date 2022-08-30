import os

from .common import *  # noqa

exacute_backend = os.environ.get("EXAQUTE_BACKEND")
print("EXAQUTE_BACKEND=", exacute_backend)

if exacute_backend:
    exacute_backend = exacute_backend.lower()

if not exacute_backend:
    print("ExaQUte backend: Local")
    from .local import *  # noqa
elif exacute_backend == "hyperloom":
    from .hyperloom import *  # noqa
elif exacute_backend == "pycompss":
    from .pycompss import *  # noqa
else:
    raise Exception("Unknown exaqute backend: {}".format(exacute_backend))
