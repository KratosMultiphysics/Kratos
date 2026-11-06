"""Processes that write solver data out as machine-learning training data.

============================== ================================================
Module                         Writes
============================== ================================================
``dataset_export_process``     ``.npz`` samples of nodal/elemental/Gauss fields
``grid_dataset_export_process``the same, resampled onto a regular voxel grid
``cae_dataset_export_process`` per-case ``.npz`` in ``physicsnemo.datapipes.cae``
                               layout (STL surface + volume fields)
``mesh_export_process``        ``.pmsh`` mesh series for
                               ``physicsnemo.datapipes.mesh_dataset``
``curator_export_process``     AI-ready Zarr stores or VTU grids, through
                               ``physicsnemo-curator``
``usd_export_process``         one time-sampled OpenUSD stage per solve - the
                               scrubbable digital twin (``usd-core``)
``streaming_dataset_export_process``
                               nothing to disk - it feeds a live queue that
                               ``training.streaming_dataset`` trains out of
============================== ================================================

The dataset and ``DataLoader`` factories that *consume* these files live in
``training.torch_dataset``.
"""
