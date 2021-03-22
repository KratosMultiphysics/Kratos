"""Defines the framework to compute in a distributed environment.

This module provides the definition of the API used throughout the library. There are several
implementation of this API. The default one is standalone, without dependencies, and provides
no parallelisation. The other ones allow distributed computing, but depend on external
libraries. Edit this module to choose a different implementation.

"""

# Without parallisation (built in XMC)
from .localEnvironment import *

# Parallelisation, dependent on external libraries:

# PyCOMPSs API
# from exaqute.ExaquteTaskPyCOMPSs import *

# Hyperloom API
# from exaqute.ExaquteTaskHyperLoom import *
