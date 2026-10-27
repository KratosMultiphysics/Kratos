"""Kratos ``Process`` implementations: everything you attach to a solve.

Every module in this package and its subpackages exposes a ``Factory`` and is
therefore usable from ``ProjectParameters.json``:

    {
        "python_module" : "inference_process",
        "kratos_module"  : "KratosMultiphysics.PhysicsNeMoApplication.processes.inference",
        "Parameters"     : { }
    }

Subpackages
-----------
``inference``
    Deploy a trained model inside the solution loop.
``export``
    Write solver data out as training data.

Three processes sit at this level because they are neither: ``adaptive_remesh_process``
changes the mesh, ``validation_metrics_process`` measures the result, and
``adjoint_sensitivity_process`` puts an exact dJ/dX field on the model part.

Anything without a ``Factory`` belongs somewhere else - ``bridges``, ``training``,
``physics``, ``deployment`` or ``distributed``.
"""
