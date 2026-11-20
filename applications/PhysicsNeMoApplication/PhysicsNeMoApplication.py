# Application dependent names and paths
# NOTE: this module must never import torch or physicsnemo (directly or
# transitively) — importing KratosMultiphysics.PhysicsNeMoApplication must
# succeed in environments without any ML package installed. All torch and
# physicsnemo imports are performed lazily inside the specific submodules
# that need them.
from KratosMultiphysics import _ImportApplication
from KratosPhysicsNeMoApplication import *
application = KratosPhysicsNeMoApplication()
application_name = "KratosPhysicsNeMoApplication"

_ImportApplication(application, application_name)
