# MVQNConvergenceAccelerator has been renamed to IQNIMVJExplicitConvergenceAccelerator. This script is for backwards compatibility.

import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from .iqnimvj_explicit import IQNIMVJExplicitConvergenceAccelerator

def Create(settings):
    cs_tools.SettingsTypeCheck(settings)
    return IQNIMVJExplicitConvergenceAccelerator(settings)