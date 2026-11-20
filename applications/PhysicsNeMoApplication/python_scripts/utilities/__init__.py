"""Small helpers with no better home.

``tensor_adaptor_dataset_utils``
    Shared gather/scatter helpers over ``Kratos.TensorAdaptors``, used by the
    bridges and by the MPI export paths.
``array_backend_utils``
    Opt-in CuPy acceleration for the array-heavy bridge paths, with numpy as
    the default and the fallback. Also the zero-copy CuPy-to-torch handoff.
``nvtx_utils``
    Opt-in NVTX ranges around deployment hot paths, for Nsight Systems.
``shallow_water_reference``
    A numpy-only linear shallow-water integrator generating trajectories for
    the GraphCast grid recipe.
``lennard_jones_reference``
    A numpy-only velocity-Verlet molecular-dynamics integrator, the data
    source of the Lennard-Jones particle recipe when neither OpenMM nor
    ``DEMApplication`` is at hand.
``point_subsampling``
    Token budgets for point-cloud models - farthest-point, uniform and
    bounding-box selection, with the nearest-selected fill that still puts a
    prediction on every node.

If something here grows a theme, give it a package.
"""
