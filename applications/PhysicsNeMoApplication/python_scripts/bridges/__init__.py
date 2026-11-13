"""Data conversion between Kratos containers and PhysicsNeMo representations.

A bridge is stateless and knows nothing about the solution loop: it turns a
``ModelPart`` into something PhysicsNeMo can consume, and scatters a prediction
back. The ``processes`` package is what calls them at the right moment.

===================== =========================================================
Module                Kratos side  ->  PhysicsNeMo side
===================== =========================================================
``torch_bridge``      nodal/elemental/Gauss data <-> ``torch.Tensor`` (zero-copy)
``mesh_bridge``       any mesh <-> ``physicsnemo.mesh.Mesh`` (subpackage; also
                      holds mesh generation, deformation, SDFs and remeshing)
``graph_bridge``      element connectivity <-> a MeshGraphNet edge graph
``grid_bridge``       unstructured fields <-> regular voxel grids
``particle_bridge``   node clouds <-> radius/kNN proximity graphs
``mapping_bridge``    non-matching transfer via MappingApplication mappers
``calculus_bridge``   tessellated mesh -> gradient/divergence/curl/Laplacian
``rom_bridge``        RomApplication POD bases <-> reduced coordinates
``adjoint_bridge``    Kratos response functions / ``SensitivityBuilder``
                      output -> row-ordered dJ/dX arrays
``cfd_bridge``        ``ModelPart`` <-> pyvista, plus ``physicsnemo-cfd`` metrics
``curator_bridge``    a running solve -> a ``physicsnemo-curator`` ETL source
``vfgn_bridge``       sintering/AM state <-> Virtual Foundry GraphNet sequences
``vtk_bridge``        ``VtkOutputProcess`` output -> upstream's ``VTKReader``
                      datasets (training-only: no provenance, no scatter-back)
``globe_bridge``      boundary sub-model-parts -> GLOBE boundary meshes plus
                      the reference lengths its kernels need
===================== =========================================================

Serialization is not conversion: ONNX and Triton export live in ``deployment``.
"""
