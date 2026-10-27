"""Kratos as the ground-truth solver in a ``physicsnemo.active_learning`` loop.

``kratos_label_strategy``
    ``CreateKratosLabelStrategy`` - the ``LabelStrategy`` implementation that
    turns a queried sample into a real Kratos solve.
``execution_backends``
    How those solves run: in-process (small problems) or as subprocesses
    (recommended - it keeps Kratos MPI ranks and ``torch.distributed`` ranks in
    separate OS processes, fans out over ``max_parallel_jobs``, and supports
    HPC job submission through an ``srun``/``sbatch --wait`` prefix).
``query_strategies``
    Which samples to label next: ensemble disagreement, predictive entropy, or
    solver-residual scoring (the real PDE residual ranking the surrogate's weak
    spots).
``metrology``
    A ``MetrologyStrategy`` backed by the validation-metrics machinery.
``sample_io``
    The ``KratosALSample`` payload and the ``.npz`` round trip a labeling case
    comes back through.
"""
