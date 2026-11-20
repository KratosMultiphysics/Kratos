"""How a queried active-learning sample becomes a Kratos solve.

``base_backend``
    ``ExecutionBackend`` - the contract: one ``KratosALSample`` in, one
    dictionary of labeled field arrays out.
``in_process_backend``
    ``InProcessBackend`` - runs the ``AnalysisStage`` in this interpreter.
    Small problems and notebooks; shares the process with torch.
``subprocess_backend``
    ``SubprocessBackend`` - the primary mode. One Kratos process per case
    from a template directory, results harvested from the ``.npz`` files the
    case's ``dataset_export_process`` writes; keeps Kratos MPI ranks and
    ``torch.distributed`` ranks in separate OS processes, fans out over
    ``max_parallel_jobs``, and becomes an HPC submission with an
    ``srun``/``sbatch --wait`` prefix in ``run_command``.
"""
