# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import sys, os

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosMultilevelMonteCarloApplication import *
application = KratosMultilevelMonteCarloApplication()
application_name = "KratosMultilevelMonteCarloApplication"

_ImportApplication(application, application_name)

# ExaQUte API
exaqute_backend = os.environ.get("EXAQUTE_BACKEND")
if exaqute_backend:
    print("EXAQUTE_BACKEND = {}".format(exaqute_backend))
    exaqute_backend = exaqute_backend.lower()
else:
    exaqute_backend = "local"
    print("Default ExaQUte backend = {}".format(exaqute_backend))

if exaqute_backend == "local":
    print("*********************************************************")
    print("** Warning: Running a mockup version of PyCOMPSs       **")
    print("** Warning: The code will NOT RUN USING PARALLEL TASKS **")
    print("*********************************************************")
    sys.path.append(os.path.join(os.path.dirname(__file__),'PyCOMPSs'))
elif exaqute_backend in ("quake", "pycompss"):
    pass
else:
    raise ValueError(
        (
            "Unknown ExaQUte backend = {}\n"
            "Supported values are 'local', 'pycompss' and 'quake'"
        ).format(exaqute_backend)
    )

# XMC
XMC_backend = os.environ.get("XMC_BACKEND")
if XMC_backend:
    print("XMC_BACKEND = {}".format(XMC_backend))
    XMC_backend = XMC_backend.lower()
else:
    XMC_backend = "local"
    print("Default XMC backend = {}".format(XMC_backend))

if XMC_backend == "external":
    pass
elif XMC_backend == "local":
    sys.path.append(os.path.join(os.path.dirname(__file__),'XMC'))
else:
    raise ValueError(
        (
            "Unknown XMC backend = {}\n"
            "Supported values are 'local' and 'external'"
        ).format(XMC_backend)
    )
