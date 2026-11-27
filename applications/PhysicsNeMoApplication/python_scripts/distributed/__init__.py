"""Running across MPI ranks.

``distributed_utils``
    Aligns ``physicsnemo.distributed.DistributedManager`` with Kratos's
    ``DataCommunicator`` (with a loud check when the two disagree), matched
    process groups and device meshes, and ``GatherModelPartToRank0`` - the
    primitive the MPI-aware export processes reconstruct topology with.
``graph_partition_utils``
    Halo-partitioned graphs for distributed graph training: per-rank subgraphs
    whose owned sets partition the global node set exactly and whose one-hop
    neighbourhoods match a serial run. *Data* parallelism. ``HaloExchangePlan``
    and ``ExchangeHaloFeatures`` move the halo copy through torch autograd
    instead of copying it in eagerly, so a gradient reaching a halo row is
    summed back onto the rank that owns it - padded per block, because gloo
    rejects the ragged all-to-all a halo needs.
``domain_parallel_utils``
    *Domain* parallelism over ``physicsnemo.domain_parallel.ShardTensor``: one
    Kratos field or grid split across the ranks, with the halo exchange and the
    mesh-wide reductions physicsnemo's handlers provide - asserted over
    CPU/gloo, the NCCL transport untested on the reference machine.
"""
