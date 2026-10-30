"""Kratos meshes as PhysicsNeMo meshes, and back.

``tessellation``
    Splits arbitrary Kratos geometries (hexahedra, prisms, pyramids, quads,
    higher-order elements) into the simplices ``physicsnemo.mesh`` requires -
    watertight on unstructured meshes via the smallest-node-id diagonal rule.
``curved_tessellation``
    The isoparametric mode: samples the exact quadratic geometry on a
    refinement lattice, with synthetic points interpolated on gather and
    dropped on scatter-back.
``provenance``
    The map from tessellated entities back to the original Kratos ones, so a
    prediction can be written onto the entity it came from. Build it once and
    cache it - rebuilding it every step dominated the bridge benchmark.
``domain_mesh_builder``
    ``DomainMesh`` export with named boundaries taken from sub-model-parts.
``generate``
    Meshes from implicit geometry: SDF primitives and combinators, marching
    cubes, 2D loop filling, tetrahedral filling of watertight surfaces, and
    ``PopulateModelPartFromMesh`` - generated geometry as real Kratos entities.
``nurbs_sampling``
    Exact NURBS (IGA) geometry sampled on a parametric lattice - the
    isogeometric analogue of the curved mode, with the gather through the
    geometry's own basis.
``spatial``
    Signed distance fields written into ordinary nodal variables, so every
    existing gather picks them up as features.
``deformation``
    Differentiable shape parameterizations (FFD, RBF, morph, displace) with
    mesh-quality energies, writing coordinates back through the position
    tensor adaptor - plus ``DeformMesh``, the mesh-level shrinkwrap and
    Sobolev displacement filter.
``operations``
    Curvature as a nodal feature, repair, subdivision, smoothing, moments,
    and ``ExtrudeModelPart`` - a planar case swept into a 3-D one.
``sampling``
    Fields sampled at arbitrary points over the mesh's own BVH: slower than
    the Kratos locator here, but differentiable, which the locator is not.
``adaptive_remeshing``
    Surrogate-error-driven remeshing: residual score -> equidistributed size
    field -> MeshingApplication MMG. Driven in the loop by
    ``processes.adaptive_remesh_process``.
"""
